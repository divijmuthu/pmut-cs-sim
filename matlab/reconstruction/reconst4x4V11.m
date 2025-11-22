% =========================================================================
% EXPERIMENTAL DATA RECONSTRUCTION SCRIPT v3.4 (Final w/ Coherence Metric)
%
% Description:
% This is the final, fully integrated script. It replaces the visual
% diagnostic plot with the correct quantitative mutual coherence calculation
% used in previous simulation scripts.
%
% v3.4 Improvements:
% - REMOVED the qualitative visual plot for the H-matrix check.
% - ADDED the correct quantitative `analyze_coherence` helper function.
% - The script now calculates and prints the max coherence value to the
%   console after the H-matrix is generated.
% =========================================================================
clear; clc; close all;

%% ===== 1. CONFIGURATION =====
fprintf('=== Experimental Data Reconstruction Script (v3.4) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output_real_data', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s/\n', output_folder);

% --- HARDWARE CONFIGURATION (Must match experiment) ---
params.disabled_pmuts_grid = [10, 11]; 
params.fixed_rx_indices_grid = [3, 14]; 
params.num_active_tx = 8;
non_rx_grid_pos = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16];
hardware_pins   = 26:39; 
params.GRID_TO_PIN_MAP = containers.Map(non_rx_grid_pos, hardware_pins);

% --- Core Physical and Simulation Parameters ---
params.c = 343;                   % Speed of sound in air (m/s)
params.fs = 1e6;                  % Simulation sample rate (Hz)
params.pmut_width_m = 0.002;      % Width of a single pMUT element (m)
params.kerf_m = 0.005;            % Spacing between pMUT elements (m)

% --- Imaging Grid Parameters ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;

% --- Pre-processing & Gating Parameters ---
params.filter_cutoff_hz = 70000;
params.filter_order = 4;
params.gate_start_time_s = (params.target_height_m / params.c) * 0.8;
params.gate_duration_s = 3.0e-3;

% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;
params.assembly_chunk_size = 25;
params.lambda_tv_reg = 25.0;
params.rho_admm = 6.73;

%% ===== 2. LOAD ALL DATA =====
fprintf('\n--- Step 1: Loading All Experimental Data ---\n');

% --- Load and average the measured pMUT excitation waveform ---
fprintf('  Please select the FOLDER containing the exported .mat waveform files...\n');
waveform_folder = uigetdir('', 'Select Waveform Folder');
if isequal(waveform_folder, 0), error('User canceled folder selection.'); end
WAVEFORM_CHANNEL_VARIABLE = 'C'; 
[base_waveform, base_waveform_fs] = load_and_average_waveform(waveform_folder, WAVEFORM_CHANNEL_VARIABLE);
fprintf('  Base waveform loaded successfully (Fs = %.2f MHz).\n', base_waveform_fs / 1e6);

% --- Load main experimental data from HDF5 files (one by one) ---
[h5_file, h5_path] = uigetfile('*.h5', 'Select the TARGET HDF5 data file');
if isequal(h5_file, 0), error('User canceled file selection.'); end
data_filepath = fullfile(h5_path, h5_file);
[bg_h5_file, bg_h5_path] = uigetfile('*.h5', 'Select the corresponding BACKGROUND HDF5 file');
if isequal(bg_h5_file, 0), error('Background file is required for reconstruction.'); end
background_filepath = fullfile(bg_h5_path, bg_h5_file);
fprintf('  -> Target: %s\n', data_filepath);
fprintf('  -> Background: %s\n', background_filepath);

fprintf('  Loading datasets from HDF5 files...\n');
target_ch_a_raw = h5read(data_filepath, '/echo_data_ch_A');
target_ch_b_raw = h5read(data_filepath, '/echo_data_ch_B');
background_ch_a_raw = h5read(background_filepath, '/echo_data_ch_A');
background_ch_b_raw = h5read(background_filepath, '/echo_data_ch_B');
tx_pin_profiles = h5read(data_filepath, '/tx_pin_profiles');

% --- Load the corresponding CSV signal profiles ---
fprintf('  Please select the 4 corresponding CSV profile files...\n');
[delays_file, prof_path] = uigetfile('*.csv', 'Select DELAYS_US CSV file');
profiles.delays = readmatrix(fullfile(prof_path, delays_file));
[freqs_file, ~] = uigetfile('*.csv', 'Select FREQUENCIES_HZ CSV file', prof_path);
profiles.frequencies = readmatrix(fullfile(prof_path, freqs_file));
[phases_file, ~] = uigetfile('*.csv', 'Select PHASES_RAD CSV file', prof_path);
profiles.phases = readmatrix(fullfile(prof_path, phases_file));
[apods_file, ~] = uigetfile('*.csv', 'Select APODIZATIONS CSV file', prof_path);
profiles.apodizations = readmatrix(fullfile(prof_path, apods_file));

params.num_acquisitions = size(tx_pin_profiles, 2);
fprintf('  Loaded data for %d acquisitions.\n', params.num_acquisitions);

%% ===== 3. PRE-PROCESS EXPERIMENTAL DATA =====
fprintf('\n--- Step 2: Pre-processing Experimental Data ---\n');
VOLTAGE_RANGE_MV = 500.0;
RESOLUTION_BITS = 14;
max_adc_value = 2^(RESOLUTION_BITS - 1) - 1;
target_a_mv = (double(target_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
target_b_mv = (double(target_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_a_mv = (double(background_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_b_mv = (double(background_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
try
    timebase = h5readatt(data_filepath, '/', 'timebase');
catch
    fprintf('  INFO: "timebase" attribute not in HDF5 file. Defaulting to 4.\n');
    timebase = 4;
end
PICO_CLOCK_FREQ_HZ = 62.5e6;
if timebase < 3, sample_interval_s = (2^timebase) / (PICO_CLOCK_FREQ_HZ * 2);
else, sample_interval_s = (timebase - 2) / PICO_CLOCK_FREQ_HZ; end
original_fs = 1 / sample_interval_s;
fprintf('  Applying high-pass filter...\n');
[b, a] = butter(params.filter_order, params.filter_cutoff_hz / (original_fs / 2), 'high');
filtered_target_a = filtfilt(b, a, target_a_mv);
filtered_target_b = filtfilt(b, a, target_b_mv);
filtered_background_a = filtfilt(b, a, background_a_mv);
filtered_background_b = filtfilt(b, a, background_b_mv);
fprintf('  Performing background subtraction...\n');
subtracted_a = filtered_target_a - filtered_background_a;
subtracted_b = filtered_target_b - filtered_background_b;
final_a = subtracted_a - mean(subtracted_a, 2);
final_b = subtracted_b - mean(subtracted_b, 2);
fprintf('  Averaging receiver channels...\n');
processed_data = (final_a + final_b) / 2.0;
fprintf('  Resampling data to %.1f MHz...\n', params.fs / 1e6);
[p, q] = rat(params.fs / original_fs);
processed_data_resampled = resample(processed_data', p, q)';
b_vector = processed_data_resampled(:);

%% ===== 4. GENERATE AND ANALYZE H-MATRIX =====
fprintf('\n--- Step 3: Generating Digital Twin H-Matrix ---\n');
tic;
[H_raw, imaging_grid] = generate_h_matrix(params, profiles, tx_pin_profiles, base_waveform, base_waveform_fs);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

fprintf('  Aligning dimensions of H and b...\n');
num_h_samples = size(H_raw, 1) / params.num_acquisitions;
num_b_samples_raw = length(b_vector) / params.num_acquisitions;
min_samples = floor(min(num_h_samples, num_b_samples_raw));
H_reshaped = reshape(full(H_raw), [floor(num_h_samples), size(H_raw, 2), params.num_acquisitions]);
H_aligned = H_reshaped(1:min_samples, :, :);
H_un_gated = sparse(reshape(H_aligned, [min_samples * params.num_acquisitions, size(H_raw, 2)]));
b_reshaped = reshape(b_vector, [floor(num_b_samples_raw), params.num_acquisitions]);
b_aligned = b_reshaped(1:min_samples, :);

fprintf('  Applying time gate and pruning...\n');
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
fprintf('  Pruning complete. Kept %d samples per acquisition.\n', size(b_gated_list{1}, 1));

fprintf('\n--- Step 4: Analyzing Final H-Matrix Coherence ---\n');
[max_coherence, ~] = analyze_coherence(H_final);
fprintf('  Max Coherence of the gated H-matrix: %.4f\n', max_coherence);
disp('Press any key in the command window to continue to reconstruction...');
pause;

%% ===== 5. RECONSTRUCT IMAGE =====
fprintf('\n--- Step 5: Reconstructing Image via TV-ADMM ---\n');
tic;
scene_matrix_size = size(imaging_grid.X_mesh);
reconstructed_image = run_admm_reconstruction(H_final, b_final, scene_matrix_size, params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);

%% ===== 6. VISUALIZE AND SAVE =====
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

function [max_coherence, coherence_matrix] = analyze_coherence(H)
    % Calculates the mutual coherence of the measurement matrix H.
    fprintf('  Calculating coherence...');
    tic;
    % Normalize the columns of H to have unit L2 norm
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        
        % Calculate the Gram matrix (matrix of inner products)
        coherence_matrix = abs(Hn' * Hn);
        
        % Set the diagonal to zero to ignore self-coherence
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        
        % The max coherence is the largest off-diagonal element
        max_coherence = full(max(coherence_matrix(:)));
    else
        % Handle case where H is empty or all zeros
        max_coherence = 0; 
        coherence_matrix = zeros(size(H,2));
    end
    fprintf(' Done. (%.2f seconds)\n', toc);
end

function [waveform, fs] = load_and_average_waveform(folder_path, channel_var)
    mat_files = dir(fullfile(folder_path, '*.mat'));
    if isempty(mat_files), error('No .mat files found in: %s', folder_path); end
    fprintf('  Found %d waveform segments. Averaging channel "%s"...\n', length(mat_files), channel_var);
    first_file_data = load(fullfile(folder_path, mat_files(1).name));
    if ~isfield(first_file_data, channel_var), error('Channel var "%s" not in .mat file.', channel_var); end
    if ~isfield(first_file_data, 'Tinterval'), error('"Tinterval" not found in .mat file.'); end
    fs = 1 / first_file_data.Tinterval;
    waveform_data = first_file_data.(channel_var);
    if size(waveform_data, 1) == 1, waveform_data = waveform_data'; end
    accumulator = double(waveform_data);
    for i = 2:length(mat_files)
        data = load(fullfile(folder_path, mat_files(i).name));
        current_waveform = data.(channel_var);
        if size(current_waveform, 1) == 1, current_waveform = current_waveform'; end
        accumulator = accumulator + double(current_waveform);
    end
    waveform = accumulator / length(mat_files);
    waveform = waveform - mean(waveform);
    fprintf('  Averaging complete.\n');
end

function [H, imaging_grid] = generate_h_matrix(config, profiles, tx_profiles_pins, base_waveform, base_fs)
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    num_x_grid = 4; num_y_grid = 4;
    pin_to_grid_map = containers.Map('KeyType','double','ValueType','double');
    all_grid_pos = config.GRID_TO_PIN_MAP.keys;
    all_pins = config.GRID_TO_PIN_MAP.values;
    for i = 1:length(all_pins), pin_to_grid_map(all_pins{i}) = all_grid_pos{i}; end
    x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
    y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
    rx_enabled_matrix = zeros(num_y_grid, num_y_grid);
    rx_enabled_matrix(config.fixed_rx_indices_grid) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_y_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0,0,0]);
    xdc_impulse(RxAperture, ones(1,10));
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    [p, q] = rat(fs / base_fs);
    base_waveform_resampled = resample(base_waveform, p, q);
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        active_tx_pins = tx_profiles_pins(:, r_acq);
        active_grid_indices = cell2mat(values(pin_to_grid_map, num2cell(active_tx_pins)));
        tx_enabled_matrix = zeros(num_y_grid, num_y_grid);
        tx_enabled_matrix(active_grid_indices) = 1;
        tx_enabled_matrix(config.disabled_pmuts_grid) = 0;
        TxAperture = xdc_2d_array(num_y_grid, num_y_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0,0,0]);
        delays_us = profiles.delays(r_acq, :);
        apod_weights = profiles.apodizations(r_acq, :);
        xdc_excitation(TxAperture, base_waveform_resampled');
        xdc_apodization(TxAperture, 0, apod_weights);
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        xdc_free(TxAperture);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, config.fs, config.assembly_chunk_size);
    xdc_free(RxAperture);
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

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix_size, config)
    imageResolution = scene_matrix_size;
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor;
    At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
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