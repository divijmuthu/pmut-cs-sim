% =========================================================================
% REAL DATA RECONSTRUCTION SCRIPT (v2.18 - Memory Safe Assembly)
%
% Description:
% This version fixes a critical memory allocation error in the "fast"
% H-matrix assembly. It now uses a robust two-pass method: the first pass
% calculates the exact number of non-zero elements required, and the second
% pass allocates the perfect amount of memory and builds the matrix. This
% is both fast and safe from memory overflow errors.
%
% v2.18 Improvements:
% - Replaced the faulty memory estimation with a two-pass system.
% - This prevents out-of-memory errors while retaining high performance.
% =========================================================================

clear; clc; close all;

%% ===== 1. CONFIGURATION AND FILE SELECTION =====
fprintf('=== pMUT REAL DATA RECONSTRUCTION (v2.18) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- Performance & Debugging Settings ---
params.num_acquisitions_to_use = Inf; 
params.downsample_factor = 4;

% --- Select Experimental Data Files ---
[h5_file, h5_path] = uigetfile('*.h5', 'Select the TARGET HDF5 data file (data_....h5)');
if isequal(h5_file, 0), error('User canceled file selection.'); end
data_filepath = fullfile(h5_path, h5_file);

[bg_h5_file, bg_h5_path] = uigetfile('*.h5', 'Select the corresponding BACKGROUND HDF5 file (or press Cancel)');
if ~isequal(bg_h5_file, 0)
    background_filepath = fullfile(bg_h5_path, bg_h5_file);
    fprintf('Loaded TARGET data: %s\n', data_filepath);
    fprintf('Loaded BACKGROUND data: %s\n', background_filepath);
else
    background_filepath = '';
    fprintf('Loaded TARGET data: %s\n', data_filepath);
    fprintf('No background file selected.\n');
end

% --- Automatically find matching profile files ---
fprintf('\n--- Automatically locating matching CSV profile files... ---\n');
[~, h5_name, ~] = fileparts(h5_file);
timestamp_match = regexp(h5_name, '\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', 'match');
if isempty(timestamp_match), error('Could not extract a valid timestamp from HDF5 filename.'); end
timestamp_str = timestamp_match{1};
fprintf('  Found timestamp in HDF5 filename: %s\n', timestamp_str);

freq_file = ['frequencies_' timestamp_str '.csv'];
delay_file = ['delays_' timestamp_str '.csv'];
phase_file = ['phases_' timestamp_str '.csv'];
cycles_file = ['cycles_' timestamp_str '.csv'];
freq_filepath = fullfile(h5_path, freq_file);
delay_filepath = fullfile(h5_path, delay_file);
phase_filepath = fullfile(h5_path, phase_file);
cycles_filepath = fullfile(h5_path, cycles_file);

if ~exist(freq_filepath, 'file') || ~exist(delay_filepath, 'file') || ~exist(phase_filepath, 'file') || ~exist(cycles_filepath, 'file')
    error('Could not automatically find all matching CSV files (including cycles.csv).');
end
fprintf('Loaded profiles:\n  - Frequencies: %s\n  - Delays: %s\n  - Phases: %s\n  - Cycles: %s\n', freq_file, delay_file, phase_file, cycles_file);

% --- Core Physical and Simulation Parameters ---
params.c = 343;
params.pmut_width_m = 0.020;
params.kerf_m = 0.001;
params.excitation_amplitude = 1e9;
params.filter_cutoff_hz = 70e3;
params.num_active_tx = 3;

% --- Imaging Grid Parameters (Top-Down XY Plane) ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.002;

% --- Time Gating Parameters ---
params.gate_start_time_s = (params.target_height_m / params.c) * 0.8;
params.gate_duration_s = 8e-3;

% --- ADMM Reconstruction Parameters ---
params.rho_admm = 7.0;
params.lambda_tv_reg = 1.5;
params.admm_tol = 1e-5;
params.admm_max_iter = 75;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;


%% ===== 2. LOAD & PRE-PROCESS EXPERIMENTAL DATA =====
fprintf('\n--- Step 1: Loading and Pre-processing Experimental Data ---\n');
tic;

target_data_raw = h5read(data_filepath, '/echo_data_raw_adc');
exp_attrs = h5info(data_filepath, '/echo_data_raw_adc').Attributes;
timebase_attr = exp_attrs(strcmp({exp_attrs.Name}, 'timebase'));
timebase = double(timebase_attr.Value);

if isempty(timebase) || ~isnumeric(timebase) || timebase < 0
    error('Could not read a valid "timebase" attribute from HDF5 file.');
end

PICO_CLOCK_FREQ_HZ = 62.5e6;
if timebase < 3, sample_interval_s = (2^timebase) / (PICO_CLOCK_FREQ_HZ * 2);
else, sample_interval_s = (timebase - 2) / PICO_CLOCK_FREQ_HZ; end
original_fs = 1 / sample_interval_s;
fprintf('  Detected original experimental sampling rate: %.2f MHz\n', original_fs / 1e6);

params.fs = original_fs / params.downsample_factor;
fprintf('  Effective sampling rate after downsampling by %dx: %.2f MHz\n', params.downsample_factor, params.fs / 1e6);

target_data = target_data_raw';

if ~isempty(background_filepath)
    fprintf('  Performing background subtraction...\n');
    background_data_raw = h5read(background_filepath, '/echo_data_raw_adc');
    background_data = background_data_raw';
    if ~isequal(size(target_data), size(background_data)), error('Target and Background data files have different dimensions after transpose.'); end
    processed_data = double(target_data) - double(background_data);
else
    processed_data = double(target_data);
end

fprintf('  Applying %.0f kHz low-pass filter...\n', params.filter_cutoff_hz / 1000);
try
    filtered_data = lowpass(processed_data, params.filter_cutoff_hz, original_fs);
catch ME, error('Signal Processing Toolbox license is required for the lowpass filter.'); end

if params.downsample_factor > 1
    fprintf('  Downsampling data by a factor of %d...\n', params.downsample_factor);
    filtered_data = downsample(filtered_data, params.downsample_factor);
end

all_freqs = double(readmatrix(freq_filepath, 'NumHeaderLines', 1));
all_delays = double(readmatrix(delay_filepath, 'NumHeaderLines', 1));
all_phases = double(readmatrix(phase_filepath, 'NumHeaderLines', 1));
all_cycles = double(readmatrix(cycles_filepath, 'NumHeaderLines', 1));

num_acqs_available = size(filtered_data, 2);
if params.num_acquisitions_to_use > num_acqs_available
    params.num_acquisitions = num_acqs_available;
else
    params.num_acquisitions = params.num_acquisitions_to_use;
end

fprintf('  Using %d out of %d available acquisitions.\n', params.num_acquisitions, num_acqs_available);
final_data = filtered_data(:, 1:params.num_acquisitions);
params.exp_frequencies_hz = all_freqs(1:params.num_acquisitions, :);
params.exp_delays_us = all_delays(1:params.num_acquisitions, :);
params.exp_phases_deg = all_phases(1:params.num_acquisitions, :);
params.exp_cycles = all_cycles(1:params.num_acquisitions, :);

b_vector = final_data(:);
fprintf('Data loading and pre-processing complete. Time: %.2f seconds.\n', toc);


%% ===== 3. GENERATE & ANALYZE H-MATRIX =====
fprintf('\n--- Step 2: Generating and Analyzing Digital Twin H-Matrix ---\n');
tic;
[H_raw, imaging_grid] = generate_h_matrix_from_profiles(params);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

% --- Align H and b vectors ---
num_h_samples_per_acq = size(H_raw, 1) / params.num_acquisitions;
num_b_samples_per_acq = size(final_data, 1);
if num_h_samples_per_acq ~= num_b_samples_per_acq
    fprintf('  Aligning H and b vectors (H samples: %d, b samples: %d)\n', ...
            num_h_samples_per_acq, num_b_samples_per_acq);
    min_samples = min(num_h_samples_per_acq, num_b_samples_per_acq);
    H_reshaped = reshape(full(H_raw), [num_h_samples_per_acq, size(H_raw, 2), params.num_acquisitions]);
    b_reshaped = reshape(b_vector, [num_b_samples_per_acq, 1, params.num_acquisitions]);
    H_aligned = H_reshaped(1:min_samples, :, :);
    b_aligned = b_reshaped(1:min_samples, :, :);
    H_un_gated = sparse(reshape(H_aligned, [min_samples * params.num_acquisitions, size(H_raw, 2)]));
    b_vector = b_aligned(:);
    fprintf('  Alignment complete. New sample length per acq: %d\n', min_samples);
else
    H_un_gated = H_raw;
end

% --- H-Matrix Diagnostics ---
fprintf('\n  --- H-Matrix Coherence Analysis ---\n');
[max_coherence_before, coherence_matrix_before] = analyze_coherence(H_un_gated);
fprintf('  Max Coherence BEFORE time-gating: %.4f\n', max_coherence_before);

% --- Apply Time Gating to BOTH b and H ---
fprintf('  Applying time gate from %.2f ms to %.2f ms...\n', ...
    params.gate_start_time_s * 1000, (params.gate_start_time_s + params.gate_duration_s) * 1000);
num_samples = size(H_un_gated, 1) / params.num_acquisitions;
time_axis = (0:num_samples-1)' / params.fs;
gate_mask_single = (time_axis >= params.gate_start_time_s) & ...
                   (time_axis < (params.gate_start_time_s + params.gate_duration_s));
gate_mask_full = repmat(gate_mask_single, params.num_acquisitions, 1);

b_gated = b_vector .* gate_mask_full;
H_gated = H_un_gated .* gate_mask_full;

[max_coherence_after, coherence_matrix_after] = analyze_coherence(H_gated);
fprintf('  Max Coherence AFTER time-gating: %.4f\n', max_coherence_after);

figure('Visible', 'on', 'Position', [100, 100, 1200, 500]);
subplot(1, 2, 1); imagesc(coherence_matrix_before); axis square; colorbar;
title(sprintf('Coherence Matrix BEFORE Gating\nMax: %.4f', max_coherence_before));
xlabel('Column Index'); ylabel('Column Index');
subplot(1, 2, 2); imagesc(coherence_matrix_after); axis square; colorbar;
title(sprintf('Coherence Matrix AFTER Gating\nMax: %.4f', max_coherence_after));
xlabel('Column Index');
saveas(gcf, fullfile(output_folder, 'coherence_analysis.png'));

fprintf('\n--- H-Matrix Analysis Complete ---\n');
fprintf('To run the full reconstruction, uncomment the final two sections of the script.\n');

%% ===== 4. RECONSTRUCT IMAGE FROM REAL DATA (COMMENTED OUT) =====
% fprintf('\n--- Step 3: Reconstructing Image via ADMM ---\n');
% tic;
% reconstructed_image = run_admm_reconstruction(H_gated, b_gated, zeros(size(imaging_grid.X_mesh)), params);
% fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);
% 
% 
% %% ===== 5. ANALYZE AND SAVE RESULTS (COMMENTED OUT) =====
% fprintf('\n--- Step 4: Analyzing and Saving Results ---\n');
% analyze_and_plot_reconstruction(reconstructed_image, imaging_grid, max_coherence_after, output_folder);
% fprintf('\n\n=== RECONSTRUCTION COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix_from_profiles(config)
    fs = config.fs; c = config.c; num_active_tx = config.num_active_tx;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    num_x_grid = 9; num_y_grid = 9; total_elements = num_x_grid * num_y_grid;
    center_offset_x = (num_x_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    center_offset_y = (num_y_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    element_centers = zeros(total_elements, 3);
    for iy = 1:num_y_grid, for ix = 1:num_x_grid
        element_no = (iy-1)*num_x_grid + ix;
        x_pos = (ix-1)*(config.pmut_width_m + config.kerf_m) - center_offset_x;
        y_pos = (iy-1)*(config.pmut_width_m + config.kerf_m) - center_offset_y;
        element_centers(element_no, :) = [x_pos, y_pos, 0];
    end, end
    tx_desired_positions = [ 25e-3, 0, 0; -12.5e-3, 21.651e-3, 0; -12.5e-3, -21.651e-3, 0 ];
    tx_active_indices = zeros(1, num_active_tx);
    for i = 1:num_active_tx, [~, tx_active_indices(i)] = min(sum((element_centers - tx_desired_positions(i,:)).^2, 2)); end
    tx_enabled_matrix = zeros(num_y_grid, num_x_grid); tx_enabled_matrix(tx_active_indices) = 1;
    TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
    [~, rx_active_index] = min(sum(element_centers.^2, 2));
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid); rx_enabled_matrix(rx_active_index) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
    y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
    all_h_data = cell(config.num_acquisitions, 1); all_start_times = zeros(config.num_acquisitions, 1); all_K_values = zeros(config.num_acquisitions, 1);
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        if mod(r_acq, 10) == 0 || r_acq == config.num_acquisitions, fprintf('  Generating H-matrix acquisition %d/%d (%.0f%%)...\n', r_acq, config.num_acquisitions, (r_acq/config.num_acquisitions)*100); end
        tx_frequencies = repmat(config.exp_frequencies_hz(r_acq), 1, num_active_tx);
        delays_us = config.exp_delays_us(r_acq, :);
        phases_deg = config.exp_phases_deg(r_acq, :);
        cycles = config.exp_cycles(r_acq, :);
        tx_signals = cell(1, num_active_tx); max_len = 0;
        for k = 1:num_active_tx
            duration = cycles(k) / tx_frequencies(k); t = 0:1/fs:duration;
            phase_rad = phases_deg(k) * pi / 180;
            apodization_factor = cycles(k) / max(config.exp_cycles(:));
            signal_base = sin(2 * pi * tx_frequencies(k) * t + phase_rad) * apodization_factor;
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx, sig = tx_signals{k}; composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig; end
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_excitation(TxAperture, ones(1, num_active_tx));
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r; all_start_times(r_acq) = start_time_r; all_K_values(r_acq) = size(h_r, 1);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    % --- *** NEW: Fast and Memory-Safe H-Matrix Assembly *** ---
    disp('  Assembling H-matrix using two-pass triplet method...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), H = sparse(config.num_acquisitions, N_pixels); return; end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices)); max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    
    % --- Pass 1: Count the exact number of non-zero elements ---
    total_nnz = 0;
    all_h_aligned = cell(config.num_acquisitions, 1);
    wb_count = waitbar(0, 'Pass 1: Counting non-zero elements...');
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned = interp1(t_current, all_h_data{r_acq}, t_common_axis, 'linear', 0);
            total_nnz = total_nnz + nnz(h_aligned);
            all_h_aligned{r_acq} = h_aligned; % Store for the next pass
        end
        waitbar(r_acq / config.num_acquisitions, wb_count);
    end
    close(wb_count);
    
    % --- Pass 2: Pre-allocate perfectly and fill the lists ---
    I_list = zeros(total_nnz, 1);
    J_list = zeros(total_nnz, 1);
    S_list = zeros(total_nnz, 1);
    currentIndex = 1;
    wb_fill = waitbar(0, 'Pass 2: Building sparse matrix triplets...');
    for r_acq = 1:config.num_acquisitions
        if ~isempty(all_h_aligned{r_acq})
            [row_idx, col_idx, s_vals] = find(all_h_aligned{r_acq});
            global_row_idx = row_idx + (r_acq - 1) * K_global_per_acq;
            num_elements = length(s_vals);
            I_list(currentIndex : currentIndex + num_elements - 1) = global_row_idx;
            J_list(currentIndex : currentIndex + num_elements - 1) = col_idx;
            S_list(currentIndex : currentIndex + num_elements - 1) = s_vals;
            currentIndex = currentIndex + num_elements;
        end
        waitbar(r_acq / config.num_acquisitions, wb_fill);
    end
    close(wb_fill);
    
    % Construct the sparse matrix in one efficient call
    total_rows = K_global_per_acq * config.num_acquisitions;
    H = sparse(I_list, J_list, S_list, total_rows, N_pixels);
    
    xdc_free(TxAperture); xdc_free(RxAperture); field_end();
end

function [max_coherence, coherence_matrix] = analyze_coherence(H)
    col_norms = vecnorm(H, 2, 1); non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0; coherence_matrix = zeros(size(H,2));
    end
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix_size, config)
    imageResolution = size(scene_matrix_size);
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor; At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    fprintf('  Starting ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 5) == 0 || k_admm == config.admm_max_iter, fprintf('  ADMM iteration %d/%d (%.0f%%)...\n', k_admm, config.admm_max_iter, (k_admm/config.admm_max_iter)*100); end
        x_prev = x_admm_img_iter;
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = admm_iteration(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            config.rho_admm, config.lambda_tv_reg, Hfun_pcg_admm, config);
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol, fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        end
        waitbar(k_admm / config.admm_max_iter, wb_admm, sprintf('ADMM Iteration %d/%d', k_admm, config.admm_max_iter));
    end
    close(wb_admm);
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

function [x_new, z_new, u_new] = admm_iteration(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, Hfun_pcg_admm, config)
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm = max(v_norm, eps);
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    u_new = u_old + opDx_tv(x_new) - z_new;
end

function analyze_and_plot_reconstruction(reconstructed_image, imaging_grid, max_coherence, output_folder)
    fprintf('  Final Max Coherence of Gated H-Matrix: %.4f\n', max_coherence);
    figure('Visible', 'on', 'Position', [100, 100, 700, 600]);
    imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Y_mesh(:,1)*1000, reconstructed_image);
    colormap(gray); colorbar; axis image;
    title(sprintf('Reconstructed Image from Real Data\n(Gated Coherence: %.4f)', max_coherence));
    xlabel('X Position (mm)'); ylabel('Y Position (mm)');
    saveas(gcf, fullfile(output_folder, 'final_reconstruction.png'));
    metrics = struct('max_coherence_gated', max_coherence);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
    save(fullfile(output_folder, 'reconstruction_data.mat'), 'reconstructed_image', 'imaging_grid');
end
