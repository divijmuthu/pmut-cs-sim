% =========================================================================
% UNIFIED PMUT ULTRASOUND SIMULATION SCRIPT (v3.2 - Final Production Version)
%
% Description:
% This is the final, definitive, all-in-one simulation script. It integrates
% the greedy profile optimization with the main imaging pipeline and includes
% full console feedback during execution. It generates the optimal CSV profiles
% needed to run the physical experiment and provides the ideal digital twin
% for reconstructing the results.
% =========================================================================
clear; clc; close all;
%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== pMUT SIMULATION (v3.2 - Final Production Version) ===\n');

% --- Profile Optimization Parameters ---
params.candidate_pool_size = 10000;

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
main_output_folder = fullfile('integrated_greedy_output_final', timestamp);
if ~exist(main_output_folder, 'dir'), mkdir(main_output_folder); end
fprintf('Saving all results to: %s/\n', main_output_folder);

% --- HARDWARE CONFIGURATION ---
params.disabled_pmut = 3;

% --- ADMM HYPERPARAMETER TUNING ---
params.lambda_tv_reg = 25.0;
params.rho_admm = 6.73;

% --- Core Physical and Simulation Parameters ---
params.c = 343;
params.fs = 1e6;
params.pmut_width_m = 0.002;
params.kerf_m = 0.005;
params.excitation_amplitude = 1e12;
params.target_SNR_db = 35;
% --- Acquisition Parameters ---
params.num_acquisitions = 150;
params.num_active_tx = 8;
params.num_active_rx = 2;
params.f_min_hz = 45e3;
params.f_max_hz = 65e3;
params.max_delay_rand_us = 500;
params.max_phase_rad = 2*pi;
params.min_apod = 0.1;
params.max_apod = 1.0;
params.use_phase_shifts = true;
params.assembly_chunk_size = 25;
% --- Imaging Grid Parameters ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;
% --- Target Scene Sweep Parameters ---
target_patterns_to_test = {'3x3_grid', '2x2_grid', 'line_5', 'cross_5'};
params.grid_spacing_mm = 25;
params.target_sizes_mm = [3, 3, 4, 4, 4, 5, 5, 5, 6];
params.target_intensities = [0.8, 0.9, 1.0, 1.0, 1.1, 1.1, 1.2, 1.2, 1.3];
% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;

%% ===== 2. GENERATE OPTIMAL PROFILES (GREEDY SEARCH) =====
fprintf('\n--- Generating Optimal Acquisition Profiles via Greedy Search... ---\n');
tic;
[profiles, profiles_table] = generate_greedy_profiles(params);

writetable(profiles_table.delays, fullfile(main_output_folder, 'profiles_delays_us.csv'));
writetable(profiles_table.frequencies, fullfile(main_output_folder, 'profiles_frequencies_hz.csv'));
writetable(profiles_table.apodizations, fullfile(main_output_folder, 'profiles_apodizations.csv'));
writetable(profiles_table.phases, fullfile(main_output_folder, 'profiles_phases_rad.csv'));
fprintf('Greedy profile generation complete. Time: %.2f seconds.\n', toc);
fprintf('Profiles saved to CSV files in output folder.\n');

%% ===== 3. GENERATE H-MATRIX USING OPTIMIZED PROFILES =====
fprintf('\n--- Generating Full H-Matrix using Optimized Profiles ---\n');
tic;
[H, imaging_grid] = generate_h_matrix(params, profiles);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);
fprintf('\n--- Analyzing H-Matrix Coherence ---\n');
[max_coherence, ~] = analyze_coherence(H);
fprintf('  Max Coherence: %.4f\n', max_coherence);

%% ===== 4. SEQUENTIAL SCENE SIMULATION SWEEP =====
for i = 1:length(target_patterns_to_test)
    params.target_pattern = target_patterns_to_test{i};
    fprintf('\n============================================================\n');
    fprintf('--- Running Simulation for Scene %d/%d: "%s" ---\n', i, length(target_patterns_to_test), params.target_pattern);
    fprintf('============================================================\n');
    
    run_output_folder = fullfile(main_output_folder, sprintf('run_%02d_%s', i, params.target_pattern));
    if ~exist(run_output_folder, 'dir'), mkdir(run_output_folder); end
    
    scene_matrix = create_scene(imaging_grid, params);
    b_vector = run_forward_simulation(H, scene_matrix(:), params);
    reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, params);
    analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, run_output_folder);
end
fprintf('\n\n=== SIMULATION SWEEP COMPLETE ===\n');

%% ===== HELPER FUNCTIONS =====

function [profiles, profiles_table] = generate_greedy_profiles(config)
    fprintf('  Step 1: Generating candidate pool of %d random profiles...\n', config.candidate_pool_size);
    tic;
    pool_delays = config.max_delay_rand_us * rand(config.candidate_pool_size, config.num_active_tx);
    pool_freqs = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(config.candidate_pool_size, config.num_active_tx);
    pool_phases = config.max_phase_rad * rand(config.candidate_pool_size, config.num_active_tx);
    pool_apods = config.min_apod + (config.max_apod - config.min_apod) * rand(config.candidate_pool_size, config.num_active_tx);
    fprintf('  ... done. Time: %.2f seconds.\n', toc);

    fprintf('  Step 2: Simulating acoustic "fingerprint" for each profile...\n');
    tic;
    fingerprint_point = [0, 0, config.target_height_m];
    fingerprints = simulate_fingerprints(pool_delays, pool_freqs, pool_phases, pool_apods, fingerprint_point, config);
    fprintf('  ... done. Time: %.2f seconds.\n', toc);

    fprintf('  Step 3: Performing greedy selection for %d optimal profiles...\n', config.num_acquisitions);
    tic;
    pool_indices = 1:config.candidate_pool_size;
    optimal_indices = zeros(config.num_acquisitions, 1);
    optimal_fingerprints = zeros(size(fingerprints, 1), config.num_acquisitions);

    optimal_indices(1) = pool_indices(1);
    optimal_fingerprints(:, 1) = fingerprints(:, 1);
    pool_indices(1) = [];

    wb_greedy = waitbar(0, 'Greedy Selection...');
    for i = 2:config.num_acquisitions
        candidates_fingerprints = fingerprints(:, pool_indices);
        current_optimal_set = optimal_fingerprints(:, 1:i-1);
        
        full_corr_matrix = corrcoef([candidates_fingerprints, current_optimal_set]);
        
        num_candidates = size(candidates_fingerprints, 2);
        corr_matrix = abs(full_corr_matrix(1:num_candidates, num_candidates+1:end));
        
        max_corrs_per_candidate = max(corr_matrix, [], 2);
        [~, min_idx_in_pool] = min(max_corrs_per_candidate);
        
        best_candidate_global_idx = pool_indices(min_idx_in_pool);
        
        optimal_indices(i) = best_candidate_global_idx;
        optimal_fingerprints(:, i) = fingerprints(:, best_candidate_global_idx);
        
        pool_indices(min_idx_in_pool) = [];
        
        waitbar(i / config.num_acquisitions, wb_greedy, sprintf('Selecting profile %d/%d', i, config.num_acquisitions));
    end
    close(wb_greedy);
    fprintf('  ... done. Time: %.2f seconds.\n', toc);
    
    profiles.delays = pool_delays(optimal_indices, :);
    profiles.frequencies = pool_freqs(optimal_indices, :);
    profiles.phases = pool_phases(optimal_indices, :);
    profiles.apodizations = pool_apods(optimal_indices, :);
    
    colNames = cellstr("Tx" + (1:config.num_active_tx));
    profiles_table.delays = array2table(profiles.delays, 'VariableNames', colNames);
    profiles_table.frequencies = array2table(profiles.frequencies, 'VariableNames', colNames);
    profiles_table.phases = array2table(profiles.phases, 'VariableNames', colNames);
    profiles_table.apodizations = array2table(profiles.apodizations, 'VariableNames', colNames);
end

function fingerprints = simulate_fingerprints(delays, freqs, phases, apods, point, config)
    num_profiles = size(delays, 1);
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    num_x_grid = 4; num_y_grid = 4;
    
    active_tx_indices = [1, 2, 4, 5, 12, 13, 15, 16];
    active_rx_indices = [6, 11];
    
    if any(ismember(active_tx_indices, config.disabled_pmut)) || any(ismember(active_rx_indices, config.disabled_pmut))
        error('Fixed geometry for fingerprinting conflicts with disabled PMUT.');
    end
    
    tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    tx_enabled_matrix(active_tx_indices) = 1;
    TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
    
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(active_rx_indices) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    all_h_data = cell(num_profiles, 1);
    all_start_times = zeros(num_profiles, 1);
    all_K_values = zeros(num_profiles, 1);
    wb = waitbar(0, 'Simulating fingerprints...');
    
    for i = 1:num_profiles
        tx_frequencies = freqs(i, :);
        delays_us = delays(i, :);
        tx_phases = phases(i, :);
        apod_weights = apods(i, :);
        
        max_len = 0;
        tx_signals = cell(1, config.num_active_tx);
        for k = 1:config.num_active_tx
            f_k = tx_frequencies(k);
            phase_k = tx_phases(k);
            duration = 5 / f_k;
            t_k = 0:1/fs:duration;
            signal_k = sin(2*pi*f_k*t_k + phase_k) .* tukeywin(length(t_k), 0.4)';
            tx_signals{k} = signal_k;
            if length(t_k) > max_len, max_len = length(t_k); end
        end
        composite_waveform = zeros(1, max_len);
        for k = 1:config.num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
        end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_apodization(TxAperture, 0, apod_weights);
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, point);
        all_h_data{i} = h_r;
        all_start_times(i) = start_time_r;
        all_K_values(i) = size(h_r, 1);
        waitbar(i/num_profiles, wb);
    end
    close(wb);
    
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    
    fingerprints = zeros(length(t_common_axis), num_profiles);
    for i = 1:num_profiles
        if all_K_values(i) > 0 && ~isempty(all_h_data{i})
            t_current = all_start_times(i) + (0:(all_K_values(i) - 1)) / fs;
            fingerprints(:, i) = interp1(t_current, all_h_data{i}(:,1), t_common_axis, 'linear', 0);
        end
    end
    
    xdc_free(TxAperture); xdc_free(RxAperture); field_end();
end

function [H, imaging_grid] = generate_h_matrix(config, profiles)
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    num_x_grid = 4; num_y_grid = 4;
    total_elements_phys = num_x_grid * num_y_grid;
    all_active_elements = setdiff(1:total_elements_phys, config.disabled_pmut);
    
    x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
    y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        
        all_rx_pairs = [ (1:8)', (17-(1:8))' ];
        if ~isempty(config.disabled_pmut), rows_to_remove = any(all_rx_pairs == config.disabled_pmut, 2); all_rx_pairs(rows_to_remove, :) = []; end
        chosen_pair_idx = randi(size(all_rx_pairs, 1));
        active_rx_indices = all_rx_pairs(chosen_pair_idx, :);
        
        available_for_tx = setdiff(all_active_elements, active_rx_indices);
        
        quad_defs = {[1, 2, 5, 6], [3, 4, 7, 8], [9, 10, 13, 14], [11, 12, 15, 16]};
        quadrants_available = cell(1,4);
        for q = 1:4, quadrants_available{q} = intersect(quad_defs{q}, available_for_tx); end
        
        active_tx_indices = zeros(1, config.num_active_tx);
        elements_per_quad = config.num_active_tx / 4;
        currentIndex = 1;
        for q = 1:4
            quad_elements = quadrants_available{q};
            if length(quad_elements) < elements_per_quad, error('A quadrant does not have enough available elements.'); end
            rand_indices = randperm(length(quad_elements), elements_per_quad);
            active_tx_indices(currentIndex : currentIndex + elements_per_quad - 1) = quad_elements(rand_indices);
            currentIndex = currentIndex + elements_per_quad;
        end
        
        num_active_tx = length(active_tx_indices);
        
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix(active_tx_indices) = 1;
        TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        rx_enabled_matrix(active_rx_indices) = 1;
        RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
        
        tx_frequencies = profiles.frequencies(r_acq, :);
        delays_us = profiles.delays(r_acq, :);
        apod_weights = profiles.apodizations(r_acq, :);
        tx_phases = profiles.phases(r_acq, :);
        
        max_len = 0;
        tx_signals = cell(1, num_active_tx);
        for k = 1:num_active_tx
            f_k = tx_frequencies(k);
            phase_k = tx_phases(k);
            duration = 5 / f_k;
            t_k = 0:1/fs:duration;
            signal_k = sin(2*pi*f_k*t_k + phase_k) .* tukeywin(length(t_k), 0.4)';
            tx_signals{k} = signal_k;
            if length(t_k) > max_len, max_len = length(t_k); end
        end
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
        end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_impulse(RxAperture, ones(1,10));
        xdc_apodization(TxAperture, 0, apod_weights);
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        xdc_free(TxAperture); xdc_free(RxAperture);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, config.fs, config.assembly_chunk_size);
    field_end();
end
function H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, fs, chunk_size)
    num_acquisitions = length(all_h_data);
    valid_indices = all_K_values > 0 & ~isinf(all_start_times);
    if ~any(valid_indices), H = sparse(num_acquisitions, N_pixels); return; end
    
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    
    num_chunks = ceil(num_acquisitions / chunk_size);
    H_chunks = cell(num_chunks, 1);
    
    for c_idx = 1:num_chunks
        start_acq = (c_idx - 1) * chunk_size + 1;
        end_acq = min(c_idx * chunk_size, num_acquisitions);
        num_acqs_in_chunk = end_acq - start_acq + 1;
        
        est_nnz_per_row = 50;
        total_nnz_chunk_est = K_global_per_acq * num_acqs_in_chunk * est_nnz_per_row;
        I_list = zeros(total_nnz_chunk_est, 1); J_list = zeros(total_nnz_chunk_est, 1); S_list = zeros(total_nnz_chunk_est, 1);
        currentIndex = 1;
        
        for r_acq_local = 1:num_acqs_in_chunk
            r_acq_global = start_acq + r_acq_local - 1;
            if all_K_values(r_acq_global) > 0 && ~isinf(all_start_times(r_acq_global))
                t_current = all_start_times(r_acq_global) + (0:(all_K_values(r_acq_global) - 1)) / fs;
                h_aligned = interp1(t_current, all_h_data{r_acq_global}, t_common_axis, 'linear', 0);
                
                [row_idx, col_idx, s_vals] = find(h_aligned);
                num_elements = length(s_vals);
                
                if (currentIndex + num_elements - 1) > length(I_list)
                    I_list = [I_list; zeros(total_nnz_chunk_est, 1)]; J_list = [J_list; zeros(total_nnz_chunk_est, 1)]; S_list = [S_list; zeros(total_nnz_chunk_est, 1)];
                end
                
                global_row_idx = row_idx + (r_acq_local - 1) * K_global_per_acq;
                I_list(currentIndex : currentIndex + num_elements - 1) = global_row_idx;
                J_list(currentIndex : currentIndex + num_elements - 1) = col_idx;
                S_list(currentIndex : currentIndex + num_elements - 1) = s_vals;
                currentIndex = currentIndex + num_elements;
            end
        end
        
        I_list(currentIndex:end) = []; J_list(currentIndex:end) = []; S_list(currentIndex:end) = [];
        total_rows_chunk = K_global_per_acq * num_acqs_in_chunk;
        H_chunks{c_idx} = sparse(I_list, J_list, S_list, total_rows_chunk, N_pixels);
    end
    
    H = vertcat(H_chunks{:});
end
function [max_coherence, coherence_matrix] = analyze_coherence(H)
    fprintf('  Calculating coherence...');
    tic;
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0; coherence_matrix = zeros(size(H,2));
    end
    fprintf(' Done. (%.2f seconds)\n', toc);
end
function scene_matrix = create_scene(imaging_grid, params)
    fprintf('\n--- Creating Target Scene ---\n');
    scene_matrix = zeros(size(imaging_grid.X_mesh));
    X_mm = imaging_grid.X_mesh * 1000;
    Y_mm = imaging_grid.Y_mesh * 1000;
    grid_spacing_mm = params.grid_spacing_mm;
    
    switch params.target_pattern
        case '3x3_grid'
            positions = [];
            for row = -1:1, for col = -1:1, positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm]; end, end
        case '2x2_grid'
            positions = [];
            for row = -0.5:1:0.5, for col = -0.5:1:0.5, positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm]; end, end
        case 'line_5'
            positions = [(-2:2)' * grid_spacing_mm, zeros(5,1)];
        case 'cross_5'
            positions = [0,0; -1,0; 1,0; 0,-1; 0,1] * grid_spacing_mm;
    end
    
    num_targets = size(positions, 1);
    target_sizes = repmat(params.target_sizes_mm, 1, ceil(num_targets/length(params.target_sizes_mm)));
    target_intensities = repmat(params.target_intensities, 1, ceil(num_targets/length(params.target_intensities)));
    
    for i = 1:num_targets
        x_pos_mm = positions(i, 1); y_pos_mm = positions(i, 2);
        target_size_mm = target_sizes(i); target_intensity = target_intensities(i);
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iy_scene] = min(abs(Y_mm(:,1) - y_pos_mm));
        target_radius_pixels = round(target_size_mm / (2 * (params.grid_step_m*1000)));
        if target_radius_pixels < 1, target_radius_pixels = 1; end
        for dy = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                ix_target = ix_scene + dx; iy_target = iy_scene + dy;
                if ix_target > 0 && ix_target <= size(X_mm, 2) && iy_target > 0 && iy_target <= size(Y_mm, 1)
                    scene_matrix(iy_target, ix_target) = target_intensity;
                end
            end
        end
    end
end
function b_vector = run_forward_simulation(H, v_true_vector, params)
    fprintf('\n--- Running Forward Simulation ---\n');
    fprintf('  Performing H * v matrix multiplication (this can be slow)...');
    tic;
    Hv_signal = H * v_true_vector;
    fprintf(' Done. (%.2f seconds)\n', toc);
    signal_power = mean(Hv_signal.^2);
    noise_variance = signal_power / (10^(params.target_SNR_db / 10));
    noise = sqrt(noise_variance) * randn(size(Hv_signal));
    b_vector = Hv_signal + noise;
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config)
    fprintf('\n--- Reconstructing Image via TV-ADMM ---\n');
    imageResolution = size(scene_matrix);
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor;
    At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    fprintf('  Starting TV-ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running TV-ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 5) == 0, fprintf('  TV-ADMM iteration %d/%d...\n', k_admm, config.admm_max_iter); end
        x_prev = x_admm_img_iter;
        
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun_admm_img(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
        [x_vec_new, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter);
        x_admm_img_iter = reshape(x_vec_new, imageResolution);
        
        kap = config.lambda_tv_reg / config.rho_admm;
        v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
        v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm(v_norm < eps) = 1;
        shr = max(0, 1 - kap ./ v_norm);
        z_admm_grad_iter = v_z_upd .* shr;
        u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
        
        if k_admm > 1 && (norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps) < config.admm_tol), fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        if ishandle(wb_admm), waitbar(k_admm / config.admm_max_iter, wb_admm); end
    end
    if ishandle(wb_admm), close(wb_admm); end
    reconstructed_image = x_admm_img_iter;
end
function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    Afun_admm = @(x) A_admm * x(:); Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end
function [Dx, Dy] = difference_operators(imageSize)
   rows = imageSize(1); cols = imageSize(2); N_img_pixels = rows * cols;
   Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
   Dx( (cols-1)*rows+1 : cols*rows , :) = 0;
   Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
   Dy( rows:rows:N_img_pixels , :) = 0;
end

function metrics = analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, output_folder)
    fprintf('\n--- Analyzing and Saving Results ---\n');
    scene_norm = scene_matrix / max(scene_matrix(:) + eps);
    recon_norm = reconstructed_image / max(abs(reconstructed_image(:)) + eps);
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    [max_coherence, ~] = analyze_coherence(H);
    
    fprintf('  PSNR: %.2f dB\n', psnr);
    fprintf('  Max Coherence: %.4f\n', max_coherence);
    
    figure('Visible','off', 'Position', [100, 100, 1200, 500]);
    subplot(1,2,1); imagesc(scene_matrix); axis image; colormap gray; title('Ground Truth');
    subplot(1,2,2); imagesc(reconstructed_image); axis image; colormap gray; title(sprintf('Reconstruction (PSNR: %.2f dB)', psnr));
    saveas(gcf, fullfile(output_folder, 'reconstruction_comparison.png'));
    close(gcf);
    
    metrics.psnr = psnr;
end