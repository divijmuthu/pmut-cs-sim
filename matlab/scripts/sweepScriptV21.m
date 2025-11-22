% =========================================================================
% Final Synthesized H-Matrix Parameter Sweep Script v17.3
%
% Description:
% This is the definitive experiment to create a low-coherence sensing matrix.
% It introduces a physically realistic, randomly scattered imaging grid where
% a minimum distance between points, based on the acoustic wavelength, is
% enforced. This prevents numerical collisions and breaks the geometric
% correlations that dominated all previous results.
%
% v17.3 Fixes:
% - Corrected the physical packing impossibility by reducing the number of
%   random imaging points to a guaranteed-to-succeed value of 100.
% =========================================================================
clear; clc; close all;

%% --- 1. DEFINE PARAMETER SWEEP CONFIGURATION ---
% --- Debug Flag ---
debug_mode = true;

% --- Base Simulation Parameters ---
base_config.c = 343;
base_config.fs = 1e6;
base_config.pmut_width_m = 0.020;
base_config.grid_width_m = 0.150;
base_config.target_distance_m = 0.150;
base_config.grid_depth_range_m = 0.10;
base_config.num_acquisitions = 100;
base_config.excitation_amplitude = 1e12;
base_config.f_min_hz = 45e3;
base_config.f_max_hz = 65e3;
base_config.tx_pool_width_m = 0.200;
% --- NEW: Parameters for the random imaging grid ---
% *** CRITICAL FIX: Reduced number of points to a physically possible value ***
base_config.num_imaging_points = 100; % Reduced to prevent packing error
base_config.min_point_separation_m = 0.007; % Enforce 7mm minimum distance

% --- Sweep Parameters (Variables to test) ---
p.num_active_tx_sweep = {3, 5, 8};
p.delay_rand_sweep_us = {12, 50};
p.apodization_sweep = {'uniform', 'random'};


%% --- 2. SETUP OUTPUT & RUN SWEEP ---
disp('--- Initializing Final H-Matrix Sweep with Spatially Separated Random Grid ---');
session_timestamp = datestr(now, 'mmddyy_HHMMSS');
base_output_dir = fullfile(pwd, 'sweep_output_separated_random_grid');
session_output_folder = fullfile(base_output_dir, session_timestamp);
if ~exist(session_output_folder, 'dir'), mkdir(session_output_folder); end
fprintf('Saving all results to: %s\n', session_output_folder);

param_grid = allcomb(p.num_active_tx_sweep, p.delay_rand_sweep_us, p.apodization_sweep);

if debug_mode
    fprintf('*** DEBUG MODE ACTIVE: Running a single test case. ***\n');
    param_grid = {8, 50, 'random'};
end

summary_results = table();
num_total_runs = size(param_grid, 1);

for i = 1:num_total_runs
    config = base_config;
    config.num_active_tx = param_grid{i, 1};
    config.max_delay_rand_us = param_grid{i, 2};
    config.apodization_mode = param_grid{i, 3};
    
    run_name = sprintf('run%03d_tx%d_del%d_%s', i, ...
        config.num_active_tx, config.max_delay_rand_us, config.apodization_mode);
    
    run_output_folder = fullfile(session_output_folder, run_name);
    mkdir(run_output_folder);
    fprintf('\n--- [%d/%d] Running Test: %s ---\n', i, num_total_runs, run_name);
    
    H = generate_h_matrix_final(config);
    diag_results = run_diagnostics(run_name, H, config, run_output_folder);
    
    current_summary = struct2table(diag_results, 'AsArray', true);
    current_summary.run_name = {run_name};
    summary_results = [summary_results; current_summary];
end

summary_filepath = fullfile(session_output_folder, 'sweep_summary.csv');
writetable(summary_results, summary_filepath);
disp('--- Parameter Sweep Complete! ---');


%% --- PRIMARY HELPER FUNCTION: H-MATRIX GENERATION (FINAL) ---
function H = generate_h_matrix_final(config)
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);
    
    % --- Define a LARGE virtual grid to act as our pool of random locations ---
    vgrid_N = 10;
    vgrid_total_elements = vgrid_N * vgrid_N;
    vgrid_pitch = config.tx_pool_width_m / (vgrid_N - 1);
    vgrid_kerf = vgrid_pitch - config.pmut_width_m;
    if vgrid_kerf < 0.0001
        error('pMUT width is too large for the virtual grid.');
    end
    TxAperture = xdc_2d_array(vgrid_N, vgrid_N, config.pmut_width_m, config.pmut_width_m, vgrid_kerf, vgrid_kerf, ones(vgrid_N), 1, 1, [0 0 0]);

    % --- Define the FIXED Receiver at the center ---
    RxAperture = xdc_2d_array(1, 1, config.pmut_width_m, config.pmut_width_m, 0, 0, ones(1,1), 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    % --- NEW: Define a RANDOMLY SCATTERED imaging grid with minimum separation ---
    fprintf('Generating random imaging grid with minimum point separation...\n');
    grid_points = zeros(config.num_imaging_points, 3);
    min_dist_sq = config.min_point_separation_m^2;
    
    max_attempts = config.num_imaging_points * 500; % Increased buffer significantly
    attempt_count = 0;
    
    for p_idx = 1:config.num_imaging_points
        is_valid_point = false;
        while ~is_valid_point
            attempt_count = attempt_count + 1;
            if attempt_count > max_attempts
                error('Could not place all imaging points. Try reducing num_imaging_points or min_point_separation_m.');
            end
            
            % Generate a candidate point
            cand_x = (rand() - 0.5) * config.grid_width_m;
            cand_z = config.target_distance_m + (rand() - 0.5) * config.grid_depth_range_m;
            candidate_point = [cand_x, 0, cand_z];
            
            % Check distance to all previously added points
            is_too_close = false;
            if p_idx > 1
                for prev_idx = 1:p_idx-1
                    dist_sq = sum((candidate_point - grid_points(prev_idx,:)).^2);
                    if dist_sq < min_dist_sq
                        is_too_close = true;
                        break;
                    end
                end
            end
            
            if ~is_too_close
                is_valid_point = true;
                grid_points(p_idx, :) = candidate_point;
            end
        end
    end
    N_pixels = size(grid_points, 1);
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k);
            t = 0:1/fs:duration;
            random_phase = 2 * pi * rand();
            signal_base = sin(2 * pi * tx_frequencies(k) * t + random_phase);
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end

        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
        end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        
        full_apod_vector = zeros(1, vgrid_total_elements);
        full_excitation_vector = zeros(1, vgrid_total_elements);
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        full_apod_vector(active_indices) = apod_weights;
        full_excitation_vector(active_indices) = 1;
        xdc_apodization(TxAperture, 0, full_apod_vector);
        xdc_excitation(TxAperture, full_excitation_vector);
        
        full_delay_vector = zeros(1, vgrid_total_elements);
        delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        full_delay_vector(active_indices) = delays_us * 1e-6;
        xdc_focus_times(TxAperture, 0, full_delay_vector);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    disp('Assembling H-matrix using interpolation...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(config.num_acquisitions, N_pixels); return;
    end
    
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    
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
    
    xdc_free(TxAperture);
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
    
    fig_coh_mat = figure('Name', ['Coherence Matrix: ' name], 'visible', 'off');
    if ~isempty(coherence_matrix)
        imagesc(coherence_matrix); colorbar;
        title(sprintf('Coherence Matrix (Max: %.2e)', max_coherence));
        xlabel('Column Index'); ylabel('Column Index');
    end
    saveas(fig_coh_mat, fullfile(output_folder, ['diag_' lower(name) '_coherence_matrix.png']));
    close(fig_coh_mat);
    
    fig_coh_hist = figure('Name', ['Coherence Hist: ' name], 'visible', 'off');
    if ~isempty(coherence_matrix)
        histogram(coherence_matrix(:), 100, 'Normalization', 'probability');
    end
    title(sprintf('Histogram of Mutual Coherence (Max: %.2e)', max_coherence));
    xlabel('Coherence Value'); ylabel('Probability'); grid on;
    saveas(fig_coh_hist, fullfile(output_folder, ['diag_' lower(name) '_coherence_histogram.png']));
    close(fig_coh_hist);
    
    fclose(logfid);
    
    results.max_coherence = max_coherence;
    results.mean_coherence = mean_coherence;
    results.sparsity = sparsity;
end
