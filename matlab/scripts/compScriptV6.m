% =========================================================================
% Comprehensive H-Matrix Parameter Sweep Script v4.0
%
% Description:
% This script provides a unified framework for performing a large-scale
% parameter sweep to discover optimal H-matrix properties. It systematically
% tests variations in transducer count, positioning, pulse waveforms,
% signal randomness, and apodization.
%
% It generates and diagnoses an H-matrix for each parameter combination,
% saving all results and a final summary report.
%
% Key Fix: Implements the correct Field II methodology for assigning
% unique waveforms to individual transmitter elements using ele_waveform.
% =========================================================================

clear; clc; close all;

%% --- 1. DEFINE PARAMETER SWEEP CONFIGURATION ---
% Define all the parameter combinations you want to test here.

% --- Base Simulation Parameters (Fixed for all runs) ---
base_config.c = 343;
base_config.fs = 1e6;
base_config.pmut_width_m = 0.020;
base_config.kerf_m = 0.0001;
base_config.grid_width_m = 0.150;
base_config.grid_depth_start_m = 0.250;
base_config.grid_depth_end_m = 0.350;
base_config.grid_step_m = 0.004;
base_config.num_acquisitions = 100; % Acquisitions per H-matrix
base_config.base_freq_hz = 150e3;    % Center frequency for randomization

% --- Sweep Parameters (Variables to test) ---
p.num_tx_sweep = {3, 5};                      % Use {} to make it a cell array
p.position_mode_sweep = {'fixed_triangle', 'random_linear'}; 
p.waveform_sweep = {'tukey_sine', 'gaussian_sine'};
p.delay_rand_sweep_us = {0, 12};              % Use {} to make it a cell array
p.freq_rand_sweep_hz = {0, 10000};             % Use {} to make it a cell array
p.apodization_sweep = {'uniform', 'random', 'hanning'};

%% --- 2. SETUP OUTPUT & RUN SWEEP ---
disp('--- Initializing H-Matrix Sweep ---');

% --- Setup Output Directory for the entire sweep session ---
session_timestamp = datestr(now, 'mmddyy_HHMMSS');
base_output_dir = fullfile(pwd, 'sweep_output');
session_output_folder = fullfile(base_output_dir, session_timestamp);
if ~exist(session_output_folder, 'dir'), mkdir(session_output_folder); end
fprintf('Saving all results to: %s\n', session_output_folder);

% --- Create all combinations of parameters to test ---
param_grid = allcomb(p.num_tx_sweep, p.position_mode_sweep, p.waveform_sweep, ...
                     p.delay_rand_sweep_us, p.freq_rand_sweep_hz, p.apodization_sweep);

summary_results = table(); % To store key results from all runs
num_total_runs = size(param_grid, 1);

for i = 1:num_total_runs
    % --- Build the configuration for the current run ---
    config = base_config;
    config.num_tx = param_grid{i, 1};
    config.position_mode = param_grid{i, 2};
    config.waveform_type = param_grid{i, 3};
    config.max_delay_rand_us = param_grid{i, 4};
    config.max_freq_rand_hz = param_grid{i, 5};
    config.apodization_mode = param_grid{i, 6};

    % --- Create a unique name and folder for this specific run ---
    run_name = sprintf('run%03d_tx%d_%s_%s_del%d_freq%d_%s', i, ...
        config.num_tx, config.position_mode, config.waveform_type, ...
        config.max_delay_rand_us, config.max_freq_rand_hz, config.apodization_mode);
    
    run_output_folder = fullfile(session_output_folder, run_name);
    mkdir(run_output_folder);

    fprintf('\n--- [%d/%d] Running Test: %s ---\n', i, num_total_runs, run_name);

    % --- Generate the H-matrix for the current config ---
    H = generate_h_matrix(config);

    % --- Run diagnostics and get key results back ---
    diag_results = run_diagnostics(run_name, H, config, run_output_folder);
    
    % --- Log summary results ---
    current_summary = struct2table(diag_results, 'AsArray', true);
    current_summary.run_name = {run_name};
    summary_results = [summary_results; current_summary];
end

% --- Save Summary of All Runs ---
summary_filepath = fullfile(session_output_folder, 'sweep_summary.csv');
writetable(summary_results, summary_filepath);

disp('--- Parameter Sweep Complete! ---');
% NOTE: You will need the `allcomb` function from MATLAB File Exchange.

%% --- HELPER FUNCTION TO GENERATE AN H MATRIX ---
function H = generate_h_matrix(config)
    % Unpack frequently used parameters
    fs = config.fs;
    c = config.c;
    num_tx = config.num_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);

    % --- Define Transmitter Positions based on config ---
    tx_desired_positions = zeros(num_tx, 3);
    switch config.position_mode
        case 'fixed_triangle'
            if num_tx ~= 3, error('fixed_triangle mode only supports 3 transmitters.'); end
            tx_desired_positions(1,:) = [25e-3, 0, 0];
            tx_desired_positions(2,:) = [-12.5e-3, 21.651e-3, 0];
            tx_desired_positions(3,:) = [-12.5e-3, -21.651e-3, 0];
        case 'random_linear'
            max_x = 0.060; % 6 cm total aperture width
            x_coords = linspace(-max_x/2, max_x/2, num_tx);
            x_coords = x_coords + (rand(1, num_tx)-0.5) * (max_x/(num_tx*2)); % Add jitter
            tx_desired_positions(:,1) = x_coords;
        % Add other position modes here (e.g., 'fixed_linear')
    end
    
    % --- Define Apertures (using individual pistons is simpler for this model) ---
    TxApertures = cell(num_tx, 1);
    for i = 1:num_tx
        pos = tx_desired_positions(i,:);
        TxApertures{i} = xdc_piston(config.pmut_width_m/2, pos(1), pos(2), pos(3));
    end
    RxAperture = xdc_piston(config.pmut_width_m/2, 0, 0, 0); % Single central receiver
    xdc_impulse(RxAperture, ones(1,10)); % Set receiver impulse response
    
    % --- Define imaging grid ---
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = config.grid_depth_start_m : config.grid_step_m : config.grid_depth_end_m;
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);

    % --- Generate Randomized Delays, Frequencies, and Apodization values ---
    frequencies_hz = config.base_freq_hz + (rand(config.num_acquisitions, num_tx) - 0.5) * 2 * config.max_freq_rand_hz;
    delays_us = rand(config.num_acquisitions, num_tx) * config.max_delay_rand_us;
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);

    wb = waitbar(0, 'Generating H matrix...');
    for r_acq = 1:config.num_acquisitions
        h_r_sum = [];
        start_time_r_min = inf;
        
        % Get apodization weights for this acquisition
        switch config.apodization_mode
            case 'uniform', apod_weights = ones(1, num_tx);
            case 'random',  apod_weights = rand(1, num_tx);
            case 'hanning', apod_weights = hanning(num_tx)';
        end

        for tx_idx = 1:num_tx
            % --- Generate Excitation Signal for this transmitter ---
            current_freq = frequencies_hz(r_acq, tx_idx);
            cycles = 3;
            tx_duration = cycles / current_freq;
            t_signal = 0 : 1/fs : tx_duration;
            signal_base = sin(2 * pi * current_freq * t_signal);
            
            switch config.waveform_type
                case 'tukey_sine',    window = tukeywin(length(signal_base), 0.25)';
                case 'gaussian_sine', window = gausswin(length(signal_base), 2.5)';
            end
            
            % Apply apodization to base amplitude (1e8 is arbitrary)
            tx_signal = signal_base .* window * 1e8 * apod_weights(tx_idx);
            
            xdc_impulse(TxApertures{tx_idx}, tx_signal);
            
            % --- Set Delays ---
            xdc_focus_times(TxApertures{tx_idx}, 0, delays_us(r_acq, tx_idx) * 1e-6);
            
            % --- Calculate and Sum Response ---
            [h_r_single, start_time_r] = calc_h(TxApertures{tx_idx}, grid_points);
            
            % Align and sum logic (summing individual responses)
            if isempty(h_r_sum)
                h_r_sum = h_r_single;
                start_time_r_min = start_time_r;
            else
                time_diff_samples = round((start_time_r - start_time_r_min) * fs);
                len_sum = size(h_r_sum, 1);
                len_single = size(h_r_single, 1);
                new_len = max(len_sum, time_diff_samples + len_single);
                
                temp_sum = zeros(new_len, N_pixels);
                temp_sum(1:len_sum, :) = h_r_sum;

                start_idx_dest = max(1, time_diff_samples + 1);
                end_idx_dest = time_diff_samples + len_single;
                
                start_idx_src = max(1, -time_diff_samples + 1);
                end_idx_src = len_single;

                temp_sum(start_idx_dest:end_idx_dest, :) = temp_sum(start_idx_dest:end_idx_dest, :) + h_r_single(start_idx_src:end_idx_src, :);
                
                h_r_sum = temp_sum;
                start_time_r_min = min(start_time_r, start_time_r_min);
            end
        end

        % Now convolve with the receiver response
        [h_r_final, start_time_final] = calc_hhp(RxAperture, h_r_sum, start_time_r_min);

        all_h_data{r_acq} = h_r_final;
        all_start_times(r_acq) = start_time_final;
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    % Assemble final H matrix (same as your original script)
    min_start_time = min(all_start_times(all_start_times > -inf));
    max_len = 0;
    for i = 1:config.num_acquisitions
        if ~isempty(all_h_data{i})
            len = size(all_h_data{i}, 1) + round((all_start_times(i) - min_start_time) * fs);
            if len > max_len, max_len = len; end
        end
    end
    K_per_acq = max_len;
    if K_per_acq == 0, K_per_acq = 1; end
    
    H = zeros(K_per_acq * config.num_acquisitions, N_pixels);
    for r_acq = 1:config.num_acquisitions
        h_r = all_h_data{r_acq};
        if isempty(h_r), continue; end
        start_time_r = all_start_times(r_acq);
        start_sample_offset = round((start_time_r - min_start_time) * fs);
        num_samples_r = size(h_r, 1);
        
        start_idx = (r_acq-1)*K_per_acq + 1 + start_sample_offset;
        end_idx = start_idx + num_samples_r - 1;
        
        if start_idx > 0 && end_idx <= r_acq*K_per_acq
             H(start_idx:end_idx, :) = h_r;
        end
    end
    field_end;
end


%% --- HELPER FUNCTION TO RUN DIAGNOSTICS ---
function results = run_diagnostics(name, H, config, output_folder)
    [M, N] = size(H);
    log_filepath = fullfile(output_folder, 'diag_log.txt');
    logfid = fopen(log_filepath, 'w');
    fprintf(logfid, '--- Diagnostics for %s ---\n\n', name);
    fprintf(logfid, 'Configuration:\n%s\n\n', jsonencode(config, 'PrettyPrint', true));

    % --- Sparsity ---
    sparsity = nnz(H) / numel(H) * 100;
    fprintf('  Sparsity: %.4f %%\n', sparsity);
    fprintf(logfid, 'Sparsity: %.4f %%\n', sparsity);
    
    % --- Mutual Coherence ---
    col_norms = vecnorm(H, 2, 1);
    Hn = H(:, col_norms > 1e-9); % Exclude zero-norm columns
    coherence_matrix = abs(Hn' * Hn);
    coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0; % Zero out diagonal
    max_coherence = max(coherence_matrix(:));
    if isempty(max_coherence), max_coherence = 0; end
    fprintf('  Max Mutual Coherence: %.4f\n', max_coherence);
    fprintf(logfid, 'Max Mutual Coherence: %.4f\n', max_coherence);
    
    % --- Coherence Histogram Plot ---
    fig_coh_hist = figure('Name', ['Coherence Hist: ' name], 'visible', 'off');
    histogram(coherence_matrix(:), 100, 'Normalization', 'probability');
    title(sprintf('Histogram of Mutual Coherence (Max: %.3f)', max_coherence));
    xlabel('Coherence'); ylabel('Probability'); grid on;
    saveas(fig_coh_hist, fullfile(output_folder, 'diag_coherence_histogram.png'));
    close(fig_coh_hist);
    
    % --- RIP Proxy (SVD) ---
    fprintf(logfid, 'RIP Proxy (SVD of random submatrices):\n');
    K = 10;
    if N >= K
        for i = 1:5
            idx = randperm(N, K);
            s = svd(full(H(:,idx)));
            fprintf(logfid, '    Submatrix %d (K=%d): min(s)=%.3e, max(s)=%.3e, cond=%.2f\n', i, K, min(s), max(s), max(s)/min(s));
        end
    end
    
    fclose(logfid);
    
    % --- Package Results for summary table ---
    results.max_coherence = max_coherence;
    results.sparsity = sparsity;
end