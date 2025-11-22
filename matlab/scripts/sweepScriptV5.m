% =========================================================================
% Comprehensive H-Matrix Parameter Sweep Script v6.5 (Definitive Fix)
%
% Description:
% This script provides a unified framework for performing a large-scale
% parameter sweep to discover optimal H-matrix properties.
%
% v6.5 Fixes & Features:
% - Replaced the fully random element selection with the robust "clustered"
%   approach from the user's v3.3 reference script. The script now finds
%   the grid elements closest to a fixed triangular pattern.
% - This change in physics will create a more focused, powerful acoustic
%   beam, dramatically increasing signal strength at the imaging grid.
% - This should be the definitive fix for the low column norms and the
%   resulting "max_mutual_coherence = 0" issue.
% =========================================================================

clear; clc; close all;

%% --- 1. DEFINE PARAMETER SWEEP CONFIGURATION ---

% --- Debug Flag ---
% Set to true to run a single, simplified test case for fast debugging.
% Set to false to run the full parameter sweep.
debug_mode = true;

% --- Base Simulation Parameters (Aligned with baseline) ---
base_config.c = 343;                % Speed of sound in air (m/s)
base_config.fs = 1e6;               % Sampling frequency (Hz)
base_config.pmut_width_m = 0.020;   % pMUT width (20mm)
base_config.kerf_m = 0.0001;          % Kerf (0.1mm)
base_config.grid_width_m = 0.150;   % Width of the imaging slice
base_config.target_distance_m = 0.250; % Centered at 250mm
base_config.grid_depth_range_m = 0.10; % Look +/- 5cm around target
base_config.grid_step_m = 0.004;
base_config.num_acquisitions = 100; % Reduced for faster debug
base_config.f_center_hz = 55e3;     % Center frequency for randomization (45-65kHz range)
base_config.excitation_amplitude = 1e8; % As per reference script

% --- Sweep Parameters (Variables to test) ---
% Note: num_active_tx is now fixed to 3 by the new logic, so it's removed from the sweep.
p.waveform_sweep = {'tukey_sine', 'gaussian_sine'};
p.delay_rand_sweep_us = {12, 50};
p.freq_rand_sweep_hz = {0, 10000};
p.apodization_sweep = {'uniform', 'random'};

%% --- 2. SETUP OUTPUT & RUN SWEEP ---
disp('--- Initializing H-Matrix Sweep ---');

% --- Setup Output Directory for the entire sweep session ---
session_timestamp = datestr(now, 'mmddyy_HHMMSS');
base_output_dir = fullfile(pwd, 'sweep_output');
session_output_folder = fullfile(base_output_dir, session_timestamp);
if ~exist(session_output_folder, 'dir'), mkdir(session_output_folder); end
fprintf('Saving all results to: %s\n', session_output_folder);

% --- Create all combinations of parameters to test ---
param_grid = allcomb(p.waveform_sweep, p.delay_rand_sweep_us, ...
                     p.freq_rand_sweep_hz, p.apodization_sweep);

% --- If in debug mode, override the grid with a single test case ---
if debug_mode
    fprintf('*** DEBUG MODE ACTIVE: Running a single test case emulating reference script. ***\n');
    param_grid = {'tukey_sine', 12, 10000, 'random'};
end

summary_results = table(); % To store key results from all runs
num_total_runs = size(param_grid, 1);

for i = 1:num_total_runs
    % --- Build the configuration for the current run ---
    config = base_config;
    config.waveform_type = param_grid{i, 1};
    config.max_delay_rand_us = param_grid{i, 2};
    config.max_freq_rand_hz = param_grid{i, 3};
    config.apodization_mode = param_grid{i, 4};
    config.num_active_tx = 3; % Fixed by the new logic

    % --- Create a unique name and folder for this specific run ---
    run_name = sprintf('run%03d_tx%d_%s_del%d_freq%d_%s', i, ...
        config.num_active_tx, config.waveform_type, ...
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


%% --- HELPER FUNCTION TO GENERATE AN H MATRIX ---
function H = generate_h_matrix(config)
    % Unpack frequently used parameters
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);

    % --- Define a fixed 2D grid of potential element locations ---
    num_x_grid = 9;
    num_y_grid = 9;
    total_elements = num_x_grid * num_y_grid;
    
    center_offset_x = (num_x_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    center_offset_y = (num_y_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    physical_element_centers = zeros(total_elements, 3);
    for iy = 1:num_y_grid
        y_pos_el = (iy-1)*(config.pmut_width_m + config.kerf_m) - center_offset_y;
        for ix = 1:num_x_grid
            x_pos_el = (ix-1)*(config.pmut_width_m + config.kerf_m) - center_offset_x;
            element_no = (iy-1)*num_x_grid + ix;
            physical_element_centers(element_no, :) = [x_pos_el, y_pos_el, 0];
        end
    end

    % --- Define the single-element Receiver (as per reference script) ---
    rx_pos = [0, 0, 0];
    rx_distances = sum((physical_element_centers - rx_pos).^2, 2);
    [~, rx_active_index] = min(rx_distances);
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(rx_active_index) = 1;
    
    RxAperture = xdc_2d_array(num_x_grid, num_y_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    % --- Define imaging grid (aligned with baseline) ---
    grid_depth_start = config.target_distance_m - config.grid_depth_range_m/2;
    grid_depth_end = config.target_distance_m + config.grid_depth_range_m/2;
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = grid_depth_start : config.grid_step_m : grid_depth_end;
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);

    % --- Define the fixed desired transmitter positions (from v3.3) ---
    tx_desired_positions = [
        25e-3, 0, 0;
        -12.5e-3, 21.651e-3, 0;
        -12.5e-3, -21.651e-3, 0
    ];
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);

    wb = waitbar(0, 'Generating H matrix...');
    for r_acq = 1:config.num_acquisitions
        
        % --- Find grid elements closest to desired positions ---
        tx_active_indices = zeros(num_active_tx, 1);
        for i = 1:num_active_tx
            distances = sum((physical_element_centers - tx_desired_positions(i,:)).^2, 2);
            [~, min_idx] = min(distances);
            tx_active_indices(i) = min_idx;
        end
        tx_active_indices = unique(tx_active_indices); % Ensure no duplicates
        
        % --- Create a temporary TxAperture for this acquisition ---
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix(tx_active_indices) = 1;
        TxAperture = xdc_2d_array(num_x_grid, num_y_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        % --- Define Apodization, Excitation, and Delays for this acquisition ---
        num_actually_active = length(tx_active_indices);

        switch config.apodization_mode
            case 'uniform', apod_weights = ones(1, num_actually_active);
            case 'random',  apod_weights = rand(1, num_actually_active);
        end
        xdc_apodization(TxAperture, 0, apod_weights);

        current_freq = config.f_center_hz + (rand() - 0.5) * 2 * config.max_freq_rand_hz;
        cycles = 3;
        duration = cycles / current_freq;
        t = 0 : 1/fs : duration;
        signal_base = sin(2 * pi * current_freq * t);
        
        switch config.waveform_type
            case 'tukey_sine',    window = tukeywin(length(t), 0.25)';
            case 'gaussian_sine', window = gausswin(length(t), 2.5)';
        end
        waveform_data = signal_base .* window * config.excitation_amplitude;
        
        xdc_impulse(TxAperture, waveform_data);
        xdc_excitation(TxAperture, ones(1, num_actually_active));
        
        delays_us = rand(1, num_actually_active) * config.max_delay_rand_us;
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        % --- Calculate the combined response from all active elements ---
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);

        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        
        % --- Free the temporary aperture ---
        xdc_free(TxAperture);
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    % --- Assemble final H matrix ---
    min_start_time = min(all_start_times(all_start_times > -inf & ~isnan(all_start_times)));
    if isempty(min_start_time), min_start_time = 0; end
    
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
        
        start_row = (r_acq-1)*K_per_acq + 1 + start_sample_offset;
        end_row = start_row + num_samples_r - 1;
        
        if start_row > (r_acq-1)*K_per_acq && end_row <= r_acq*K_per_acq
             H(start_row:end_row, :) = h_r;
        end
    end

    xdc_free(RxAperture);
    field_end;
end


%% --- HELPER FUNCTION TO RUN DIAGNOSTICS ---
function results = run_diagnostics(name, H, config, output_folder)
    [M, N] = size(H);
    log_filepath = fullfile(output_folder, 'diag_log.txt');
    logfid = fopen(log_filepath, 'w');
    fprintf(logfid, '--- Diagnostics for %s ---\n\n', name);
    config_str = evalc('disp(config)');
    fprintf(logfid, 'Configuration:\n%s\n\n', config_str);

    % --- Sparsity ---
    sparsity = nnz(H) / numel(H) * 100;
    fprintf('  Sparsity: %.4f %%\n', sparsity);
    fprintf(logfid, 'Sparsity: %.4f %%\n', sparsity);
    
    % --- Mutual Coherence ---
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols = col_norms > 1e-9;
    Hn = H(:, non_zero_cols); 
    
    % --- Diagnostic Printout ---
    fprintf('  Coherence Check: Found %d non-zero columns (out of %d total) for calculation.\n', size(Hn, 2), N);
    fprintf(logfid, 'Coherence Check: Found %d non-zero columns (out of %d total) for calculation.\n', size(Hn, 2), N);
    fprintf(logfid, '  Column Norm Stats: Min=%.2e, Max=%.2e, Mean=%.2e\n', min(col_norms), max(col_norms), mean(col_norms));

    if isempty(Hn) || size(Hn, 2) < 2
        max_coherence = 0;
        coherence_matrix = [];
    else
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0; % Zero out diagonal
        max_coherence = max(coherence_matrix(:));
        if isempty(max_coherence), max_coherence = 0; end
    end
    fprintf('  Max Mutual Coherence: %.4f\n', max_coherence);
    fprintf(logfid, 'Max Mutual Coherence: %.4f\n', max_coherence);
    
    % --- NEW: Coherence Matrix Visualization ---
    fig_coh_mat = figure('Name', ['Coherence Matrix: ' name], 'visible', 'off');
    if ~isempty(coherence_matrix)
        imagesc(coherence_matrix);
        colorbar;
        title(sprintf('Coherence Matrix (Max: %.3f)', max_coherence));
        xlabel('Column Index'); ylabel('Column Index');
    end
    saveas(fig_coh_mat, fullfile(output_folder, ['diag_' lower(name) '_coherence_matrix.png']));
    close(fig_coh_mat);

    % --- NEW: H-Matrix Column Visualization ---
    fig_h_cols = figure('Name', ['H Columns: ' name], 'visible', 'off');
    if size(Hn, 2) > 0
        num_cols_to_plot = min(size(Hn, 2), 5);
        indices_to_plot = round(linspace(1, size(Hn, 2), num_cols_to_plot));
        plot(Hn(:, indices_to_plot));
        title('Sample Normalized H-Matrix Columns');
        xlabel('Row Index'); ylabel('Amplitude');
        legend(arrayfun(@(x) sprintf('Col %d', x), find(non_zero_cols, num_cols_to_plot), 'UniformOutput', false));
    end
    saveas(fig_h_cols, fullfile(output_folder, ['diag_' lower(name) '_H_columns.png']));
    close(fig_h_cols);

    % --- Coherence Histogram Plot ---
    fig_coh_hist = figure('Name', ['Coherence Hist: ' name], 'visible', 'off');
    if ~isempty(coherence_matrix)
        histogram(coherence_matrix(:), 100, 'Normalization', 'probability');
    end
    title(sprintf('Histogram of Mutual Coherence (Max: %.3f)', max_coherence));
    xlabel('Coherence'); ylabel('Probability'); grid on;
    saveas(fig_coh_hist, fullfile(output_folder, ['diag_' lower(name) '_coherence_histogram.png']));
    close(fig_coh_hist);
    
    % --- RIP Proxy (SVD) ---
    fprintf(logfid, 'RIP Proxy (SVD of random submatrices):\n');
    K = 10;
    if size(Hn, 2) >= K
        for i = 1:5
            idx = randperm(size(Hn, 2), K);
            s = svd(full(Hn(:,idx)));
            min_s = min(s); max_s = max(s);
            if min_s > 0
                cond_s = max_s/min_s;
            else
                cond_s = inf;
            end
            fprintf(logfid, '    Submatrix %d (K=%d): min(s)=%.3e, max(s)=%.3e, cond=%.2f\n', i, K, min_s, max_s, cond_s);
        end
    end
    
    fclose(logfid);
    
    % --- Package Results for summary table ---
    results.max_coherence = max_coherence;
    results.sparsity = sparsity;
end
