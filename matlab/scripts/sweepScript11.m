% =========================================================================
% Final Synthesized H-Matrix Parameter Sweep Script v9.0
%
% Description:
% This definitive script integrates a sophisticated waveform generation
% model into the robust sweeping framework. For each acquisition, every
% active transducer is assigned its own unique, random frequency. These
% signals are combined into a single composite waveform, maximizing the
% diversity of the transmitted acoustic field.
%
% v9.0 Final Features:
% - Implements per-element frequency diversity based on the 'quantum' script.
% - Uses a robust composite waveform method for physical accuracy.
% - Simplifies the sweep to focus on the most impactful parameters:
%   1. Number of Active Transmitters
%   2. Maximum Time Delay
%   3. Apodization Mode
% =========================================================================
clear; clc; close all;

%% --- 1. DEFINE PARAMETER SWEEP CONFIGURATION ---
% --- Debug Flag ---
% Set to true to run a single, simplified test case for fast debugging.
% Set to false to run the full parameter sweep.
debug_mode = true;

% --- Base Simulation Parameters ---
base_config.c = 343;                % Speed of sound in air (m/s)
base_config.fs = 1e6;               % Sampling frequency (Hz)
base_config.pmut_width_m = 0.020;   % pMUT width (20mm)
base_config.kerf_m = 0.0001;          % Kerf (0.1mm)
base_config.grid_width_m = 0.150;   % Width of the imaging slice
base_config.target_distance_m = 0.250; % Centered at 250mm
base_config.grid_depth_range_m = 0.10; % Look +/- 5cm around target
base_config.grid_step_m = 0.004;
base_config.num_acquisitions = 100; % Number of random acquisitions
base_config.excitation_amplitude = 1e12; % High amplitude for strong signal
base_config.f_min_hz = 45e3;        % Minimum frequency for randomization
base_config.f_max_hz = 65e3;        % Maximum frequency for randomization

% --- Sweep Parameters (Variables to test) ---
p.num_active_tx_sweep = {3, 5, 8}; % Test different numbers of active elements
p.delay_rand_sweep_us = {12, 50};  % Test different ranges of random time delays
p.apodization_sweep = {'uniform', 'random'}; % Test different amplitude weightings


%% --- 2. SETUP OUTPUT & RUN SWEEP ---
disp('--- Initializing Final Synthesized H-Matrix Sweep ---');

% --- Setup Output Directory ---
session_timestamp = datestr(now, 'mmddyy_HHMMSS');
base_output_dir = fullfile(pwd, 'sweep_output_final');
session_output_folder = fullfile(base_output_dir, session_timestamp);
if ~exist(session_output_folder, 'dir'), mkdir(session_output_folder); end
fprintf('Saving all results to: %s\n', session_output_folder);

% --- Create all combinations of parameters to test ---
param_grid = allcomb(p.num_active_tx_sweep, p.delay_rand_sweep_us, p.apodization_sweep);

% --- If in debug mode, override the grid with a single test case ---
if debug_mode
    fprintf('*** DEBUG MODE ACTIVE: Running a single test case. ***\n');
    param_grid = {3, 12, 'random'};
end

summary_results = table(); % To store key results from all runs
num_total_runs = size(param_grid, 1);

for i = 1:num_total_runs
    % --- Build the configuration for the current run ---
    config = base_config;
    config.num_active_tx = param_grid{i, 1};
    config.max_delay_rand_us = param_grid{i, 2};
    config.apodization_mode = param_grid{i, 3};
    
    % --- Create a unique name and folder for this specific run ---
    run_name = sprintf('run%03d_tx%d_del%d_%s', i, ...
        config.num_active_tx, config.max_delay_rand_us, config.apodization_mode);
    
    run_output_folder = fullfile(session_output_folder, run_name);
    mkdir(run_output_folder);
    fprintf('\n--- [%d/%d] Running Test: %s ---\n', i, num_total_runs, run_name);
    
    % --- Generate the H-matrix using the synthesized logic ---
    H = generate_h_matrix_final(config);
    
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


%% --- PRIMARY HELPER FUNCTION: H-MATRIX GENERATION (FINAL) ---
function H = generate_h_matrix_final(config)
    % Unpack frequently used parameters
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);
    
    % --- Define a fixed 2D grid of potential element locations ---
    num_x_grid = 9; num_y_grid = 9;
    total_elements = num_x_grid * num_y_grid;
    center_offset_x = (num_x_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    center_offset_y = (num_y_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    physical_element_centers = zeros(total_elements, 3);
    for iy = 1:num_y_grid, for ix = 1:num_x_grid
        element_no = (iy-1)*num_x_grid + ix;
        x_pos_el = (ix-1)*(config.pmut_width_m + config.kerf_m) - center_offset_x;
        y_pos_el = (iy-1)*(config.pmut_width_m + config.kerf_m) - center_offset_y;
        physical_element_centers(element_no, :) = [x_pos_el, y_pos_el, 0];
    end, end

    % --- Define the single-element Receiver ---
    [~, rx_active_index] = min(sum(physical_element_centers.^2, 2));
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(rx_active_index) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    % --- Define imaging grid ---
    grid_depth_start = config.target_distance_m - config.grid_depth_range_m/2;
    grid_depth_end = config.target_distance_m + config.grid_depth_range_m/2;
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = grid_depth_start : config.grid_step_m : grid_depth_end;
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    
    % --- Data storage for interpolation-based assembly ---
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        % --- Create a temporary TxAperture for this acquisition ---
        active_indices = randperm(total_elements, num_active_tx);
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix(active_indices) = 1;
        TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        % --- NEW: Generate a COMPOSITE WAVEFORM with per-element frequency diversity ---
        % 1. Generate unique random frequencies for each active element
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        
        % 2. Create the individual signals
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k);
            t = 0:1/fs:duration;
            signal_base = sin(2 * pi * tx_frequencies(k) * t);
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end

        % 3. Sum into a single composite waveform
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
        end
        
        % 4. Set the final impulse response, scaled by amplitude
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_excitation(TxAperture, ones(1, num_active_tx));
        
        % --- Apply Apodization and Delays ---
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        xdc_apodization(TxAperture, 0, apod_weights);
        
        delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        % --- Calculate and store the response ---
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        xdc_free(TxAperture);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    % --- ASSEMBLE H-MATRIX USING INTERPOLATION ---
    disp('Assembling H-matrix using interpolation...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(config.num_acquisitions, N_pixels);
        disp('Warning: No valid acquisitions found. Returning empty H matrix.');
        return;
    end
    
    all_end_times = zeros(config.num_acquisitions, 1);
    all_end_times(valid_indices) = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times(valid_indices));
    
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0, K_global_per_acq = 1; end
    
    total_rows = K_global_per_acq * config.num_acquisitions;
    estimated_nnz = round(sum(all_K_values) * N_pixels * 0.1);
    H = spalloc(total_rows, N_pixels, estimated_nnz);
    
    current_row_offset = 0;
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned = interpolate_h_response(t_current, all_h_data{r_acq}, t_common_axis);
            row_indices = current_row_offset + (1:K_global_per_acq);
            if max(row_indices) <= size(H, 1)
                H(row_indices, :) = h_aligned;
            end
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    
    xdc_free(RxAperture);
    field_end;
end

%% --- HELPER: H-RESPONSE INTERPOLATION ---
function h_aligned = interpolate_h_response(t_current, h_current, t_common)
    N_pixels = size(h_current, 2);
    h_aligned = zeros(length(t_common), N_pixels);
    if length(t_current) > 1
        for px_col = 1:N_pixels
            h_aligned(:, px_col) = interp1(t_current, h_current(:, px_col), t_common, 'linear', 0);
        end
    end
end

%% --- HELPER FUNCTION TO RUN DIAGNOSTICS ---
function results = run_diagnostics(name, H, config, output_folder)
    [~, N] = size(H);
    log_filepath = fullfile(output_folder, 'diag_log.txt');
    logfid = fopen(log_filepath, 'w');
    fprintf(logfid, '--- Diagnostics for %s ---\n\n', name);
    config_str = evalc('disp(config)');
    fprintf(logfid, 'Configuration:\n%s\n\n', config_str);
    
    % Sparsity and H-Value Distribution
    non_zero_elements = nonzeros(H);
    sparsity = numel(non_zero_elements) / numel(H) * 100;
    fprintf('  Sparsity: %.4f %%\n', sparsity);
    fprintf(logfid, 'Sparsity: %.4f %%\n', sparsity);
    
    fig_h_hist = figure('Name', ['H Values Hist: ' name], 'visible', 'off');
    if ~isempty(non_zero_elements), histogram(non_zero_elements, 100); end
    title('Histogram of Non-Zero H-Matrix Values');
    xlabel('Amplitude'); ylabel('Frequency'); grid on;
    saveas(fig_h_hist, fullfile(output_folder, ['diag_' lower(name) '_H_values_histogram.png']));
    close(fig_h_hist);
    
    % Mutual Coherence
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    
    fprintf('  Coherence Check: Found %d non-zero columns (out of %d total).\n', size(Hn, 2), N);
    fprintf(logfid, 'Coherence Check: Found %d non-zero columns (out of %d total).\n', size(Hn, 2), N);
    
    coherence_matrix = [];
    if isempty(Hn) || size(Hn, 2) < 2
        max_coherence = 0; mean_coherence = 0;
    else
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        
        max_coherence = full(max(coherence_matrix(:)));
        mean_coherence = full(mean(coherence_matrix(:)));
        if isempty(max_coherence), max_coherence = 0; end
        if isempty(mean_coherence), mean_coherence = 0; end
    end
    fprintf('  Max Mutual Coherence: %.4e\n', max_coherence);
    fprintf(logfid, 'Max Mutual Coherence: %.4e\n', max_coherence);
    
    % Coherence Matrix Visualization
    fig_coh_mat = figure('Name', ['Coherence Matrix: ' name], 'visible', 'off');
    if ~isempty(coherence_matrix)
        imagesc(coherence_matrix); colorbar;
        title(sprintf('Coherence Matrix (Max: %.2e)', max_coherence));
        xlabel('Column Index'); ylabel('Column Index');
    end
    saveas(fig_coh_mat, fullfile(output_folder, ['diag_' lower(name) '_coherence_matrix.png']));
    close(fig_coh_mat);
    
    % Coherence Histogram Plot
    fig_coh_hist = figure('Name', ['Coherence Hist: ' name], 'visible', 'off');
    if ~isempty(coherence_matrix)
        histogram(coherence_matrix(:), 100, 'Normalization', 'probability');
    end
    title(sprintf('Histogram of Mutual Coherence (Max: %.2e)', max_coherence));
    xlabel('Coherence Value'); ylabel('Probability'); grid on;
    saveas(fig_coh_hist, fullfile(output_folder, ['diag_' lower(name) '_coherence_histogram.png']));
    close(fig_coh_hist);
    
    fclose(logfid);
    
    % Package Results for summary table
    results.max_coherence = max_coherence;
    results.mean_coherence = mean_coherence;
    results.sparsity = sparsity;
end