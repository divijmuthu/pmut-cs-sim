% =========================================================================
% EXPERIMENTAL DATA RECONSTRUCTION SCRIPT (v5.5 - Complete)
%
% Description:
% This is the complete and runnable script. It contains the correct
% 1-16 grid pin mapping and data loading procedures for your experiment.
% =========================================================================
clear; clc; close all;

%% ===== 1. CONFIGURATION =====
fprintf('=== Experimental Data Reconstruction Script (v5.5) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output_real_data', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s/\n', output_folder);

% --- HARDWARE CONFIGURATION ---
% The 12 available transmitter pins from the Python script (1-16 grid)
params.CH_PINS = [1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 15, 16];
% This map translates the pin numbers from the HDF5 file to their
% corresponding 1-16 grid positions for the simulation.
params.PIN_TO_GRID_MAP = containers.Map(params.CH_PINS, params.CH_PINS);
params.receivers_grid = [3, 14];
params.disabled_grid = [10, 11];

% --- Core Physical and Simulation Parameters ---
params.c = 343;
params.fs = 2e6;
params.pmut_width_m = 0.002;
params.kerf_m = 0.005;
params.excitation_amplitude = 1.0;

% --- Imaging Grid & Pre-processing ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;
params.filter_cutoff_hz = 70000;
params.filter_order = 4;
params.gate_start_time_s = (params.target_height_m / params.c) * 0.8;
params.gate_duration_s = 3.0e-3;

% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;
params.lambda_tv_reg = 25.0;
params.rho_admm = 6.73;

%% ===== 2. LOAD EXPERIMENTAL DATA AND PROFILES =====
fprintf('\n--- Step 1: Loading Experimental Data and Profiles ---\n');

[h5_file, h5_path] = uigetfile('*.h5', 'Select the TARGET HDF5 data file');
if isequal(h5_file, 0), error('User canceled file selection.'); end
data_filepath = fullfile(h5_path, h5_file);

[bg_h5_file, bg_h5_path] = uigetfile('*.h5', 'Select the corresponding BACKGROUND HDF5 file');
if isequal(bg_h5_file, 0), error('Background file is required.'); end
background_filepath = fullfile(bg_h5_path, bg_h5_file);

fprintf('  Loading data from HDF5 files...\n');
target_ch_a_raw = h5read(data_filepath, '/echo_data_ch_A');
target_ch_b_raw = h5read(data_filepath, '/echo_data_ch_B');
background_ch_a_raw = h5read(background_filepath, '/echo_data_ch_A');
background_ch_b_raw = h5read(background_filepath, '/echo_data_ch_B');
tx_pin_profiles = h5read(data_filepath, '/tx_pin_profiles');
params.num_acquisitions = size(tx_pin_profiles, 2);

fprintf('  Please select the 2 corresponding CSV profile files...\n');
[delays_file, prof_path] = uigetfile('*.csv', 'Select DELAYS_US CSV file');
profiles.delays = readmatrix(fullfile(prof_path, delays_file));
[freqs_file, ~] = uigetfile('*.csv', 'Select FREQUENCIES_HZ CSV file', prof_path);
profiles.frequencies = readmatrix(fullfile(prof_path, freqs_file));

fprintf('  Loaded data for %d acquisitions.\n', params.num_acquisitions);

%% ===== 3. PRE-PROCESS EXPERIMENTAL DATA =====
fprintf('\n--- Step 2: Pre-processing Experimental Data ---\n');
VOLTAGE_RANGE_MV = 5000.0;
RESOLUTION_BITS = 14;
max_adc_value = 2^(RESOLUTION_BITS - 1) - 1;
target_a_mv = (double(target_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
target_b_mv = (double(target_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_a_mv = (double(background_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_b_mv = (double(background_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;

fprintf('  Performing background subtraction...\n');
subtracted_a = target_a_mv - background_a_mv;
subtracted_b = target_b_mv - background_b_mv;

original_fs = 2.5e6;
fprintf('  Assuming original sample rate of %.2f MHz from PicoScope settings.\n', original_fs / 1e6);
fprintf('  Applying low-pass filter...\n');
[b, a] = butter(params.filter_order, params.filter_cutoff_hz / (original_fs / 2), 'low');
filtered_a = filtfilt(b, a, subtracted_a);
filtered_b = filtfilt(b, a, subtracted_b);

fprintf('  Averaging receiver channels...\n');
processed_data = (filtered_a + filtered_b) / 2.0;

fprintf('  Resampling data to %.1f MHz...\n', params.fs / 1e6);
[p, q] = rat(params.fs / original_fs);
b_vector_resampled = resample(processed_data, p, q);
b_vector = b_vector_resampled(:);

%% ===== 4. GENERATE DIGITAL TWIN H-MATRIX =====
fprintf('\n--- Step 3: Generating Digital Twin H-Matrix ---\n');
tic;
[H_raw, imaging_grid] = generate_h_matrix_digital_twin(params, profiles, tx_pin_profiles);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

% --- Align and Gate Data ---
fprintf('  Aligning and gating data...\n');
num_h_samples = size(H_raw, 1) / params.num_acquisitions;
num_b_samples_raw = length(b_vector) / params.num_acquisitions;
min_samples = floor(min(num_h_samples, num_b_samples_raw));
H_reshaped = reshape(full(H_raw), [floor(num_h_samples), size(H_raw, 2), params.num_acquisitions]);
H_aligned = H_reshaped(1:min_samples, :, :);
H_un_gated = sparse(reshape(H_aligned, [min_samples * params.num_acquisitions, size(H_raw, 2)]));
b_reshaped = reshape(b_vector, [floor(num_b_samples_raw), params.num_acquisitions]);
b_aligned = b_reshaped(1:min_samples, :);
time_axis = (0:min_samples-1)' / params.fs;
gate_mask_single = (time_axis >= params.gate_start_time_s) & (time_axis < (params.gate_start_time_s + params.gate_duration_s));
H_gated_list = cell(params.num_acquisitions, 1);
b_gated_list = cell(params.num_acquisitions, 1);
for i = 1:params.num_acquisitions
    H_single_acq = H_un_gated((i-1)*min_samples + 1 : i*min_samples, :);
    H_gated_list{i} = H_single_acq(gate_mask_single, :);
    b_gated_list{i} = b_aligned(gate_mask_single, i);
end
H_final = vertcat(H_gated_list{:});
b_final = vertcat(b_gated_list{:});

%% ===== 5. DIAGNOSTIC: COHERENCE CHECK =====
fprintf('\n--- Step 4: Analyzing Final H-Matrix Coherence ---\n');
[max_coherence, coherence_matrix] = analyze_coherence(H_final);
fprintf('  Max Coherence of the final H-matrix: %.4f\n', max_coherence);
figure('Name', 'Coherence Matrix');
imagesc(coherence_matrix);
axis image; colorbar;
title(sprintf('Coherence Matrix (Max Value: %.4f)', max_coherence));
xlabel('Pixel Index'); ylabel('Pixel Index');

%% ===== 6. RECONSTRUCT IMAGE =====
fprintf('\n--- Step 5: Reconstructing Image via TV-ADMM ---\n');
tic;
scene_matrix_size = size(imaging_grid.X_mesh);
reconstructed_image = run_admm_reconstruction(H_final, b_final, scene_matrix_size, params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);

%% ===== 7. VISUALIZE AND SAVE =====
fprintf('\n--- Step 6: Visualizing and Saving Results ---\n');
figure('Name', 'Final Reconstruction', 'Position', [100, 100, 700, 600]);
imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Y_mesh(:,1)*1000, reconstructed_image);
colormap(gray); colorbar; axis image;
title(sprintf('Final Reconstructed Image (%s)', timestamp));
xlabel('X Position (mm)'); ylabel('Y Position (mm)');
saveas(gcf, fullfile(output_folder, 'final_reconstruction.png'));
save(fullfile(output_folder, 'reconstruction_data.mat'), 'reconstructed_image', 'imaging_grid', 'params');
fprintf('Reconstruction saved to: %s\n', output_folder);

%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix_digital_twin(config, profiles, tx_pin_profiles)
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    num_x_grid = 4; num_y_grid = 4;
    
    x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
    y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
    
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(config.receivers_grid) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    BURST_CYCLES = 7;
    
    wb = waitbar(0, 'Generating Digital Twin H-Matrix...');
    for r_acq = 1:config.num_acquisitions
        active_tx_pins = tx_pin_profiles(:, r_acq);
        active_grid_indices = values(config.PIN_TO_GRID_MAP, num2cell(active_tx_pins));
        
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix([active_grid_indices{:}]) = 1;
        tx_enabled_matrix(config.disabled_grid) = 0;
        TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        tx_frequencies_hz = profiles.frequencies(r_acq, :);
        delays_us = profiles.delays(r_acq, :);
        
        schedules = cell(8, 1);
        max_total_ticks = 0;
        for i = 1:8
            delay_ticks = delays_us(i) * (fs / 1e6);
            half_period_ticks = round((fs / tx_frequencies_hz(i)) / 2);
            edges = zeros(1, BURST_CYCLES * 2);
            for j = 1:BURST_CYCLES * 2
                edges(j) = round(delay_ticks + (j-1) * half_period_ticks);
            end
            schedules{i} = edges;
            if edges(end) > max_total_ticks, max_total_ticks = edges(end); end
        end
        
        time_vector = 0:max_total_ticks;
        composite_wave = zeros(1, length(time_vector));
        for i = 1:8
            edges = schedules{i};
            level = 1;
            for j = 1:length(edges)-1
                start_idx = edges(j) + 1;
                end_idx = edges(j+1);
                if end_idx > length(composite_wave), end_idx = length(composite_wave); end
                if start_idx > end_idx, continue; end
                composite_wave(start_idx:end_idx) = composite_wave(start_idx:end_idx) + level;
                level = -level;
            end
        end
        
        xdc_excitation(TxAperture, composite_wave * config.excitation_amplitude);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        xdc_free(TxAperture);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    xdc_free(RxAperture);
    H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, config.fs, 25);
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
        H_chunk_list = cell(num_acqs_in_chunk, 1);
        for r_acq_local = 1:num_acqs_in_chunk
            r_acq_global = start_acq + r_acq_local - 1;
            if all_K_values(r_acq_global) > 0 && ~isinf(all_start_times(r_acq_global))
                t_current = all_start_times(r_acq_global) + (0:(all_K_values(r_acq_global) - 1)) / fs;
                h_aligned = interp1(t_current, all_h_data{r_acq_global}, t_common_axis, 'linear', 0);
                H_chunk_list{r_acq_local} = sparse(h_aligned);
            else
                H_chunk_list{r_acq_local} = sparse(K_global_per_acq, N_pixels);
            end
        end
        H_chunks{c_idx} = vertcat(H_chunk_list{:});
    end
    H = vertcat(H_chunks{:});
end

function [max_coherence, coherence_matrix] = analyze_coherence(H)
    fprintf('  Calculating coherence...');
    tic;
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if size(Hn, 2) > 1
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0; 
        coherence_matrix = zeros(size(H,2));
    end
    fprintf(' Done. (%.2f seconds)\n', toc);
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix_size, config)
    imageResolution = size(scene_matrix_size);
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
        if mod(k_admm, 10) == 0, fprintf('  TV-ADMM iteration %d/%d...\n', k_admm, config.admm_max_iter); end
        x_prev = x_admm_img_iter;
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun_admm_img(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
        [x_vec_new, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_prev(:));
        x_admm_img_iter = reshape(x_vec_new, imageResolution);
        kap = config.lambda_tv_reg / config.rho_admm;
        v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
        v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm(v_norm < eps) = 1;
        shr = max(0, 1 - kap ./ v_norm);
        z_admm_grad_iter = v_z_upd .* shr;
        u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
        if k_admm > 1 && (norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps) < config.admm_tol), fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        if ishandle(wb_admm), waitbar(k_admm / config.admm_max_iter, wb_admm, sprintf('ADMM Iteration %d / %d', k_admm, config.admm_max_iter)); end
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