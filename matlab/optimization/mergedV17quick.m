% =========================================================================
% REAL DATA RECONSTRUCTION SCRIPT (v2.10 - Performance Options)
%
% Description:
% This version adds performance and debugging options to manage the long
% processing times. It includes a progress bar for the slow interpolation
% step and allows for optional downsampling and using a subset of
% acquisitions for faster testing.
%
% v2.10 Improvements:
% - Added `params.num_acquisitions_to_use` to run tests on a smaller dataset.
% - Added `params.downsample_factor` to reduce data resolution for speed.
% - Added a waitbar to the H-matrix assembly loop to provide feedback.
% =========================================================================

clear; clc; close all;

%% ===== 1. CONFIGURATION AND FILE SELECTION =====
fprintf('=== pMUT REAL DATA RECONSTRUCTION (v2.10) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- Performance & Debugging Settings ---
% To speed up testing, you can use a subset of the acquisitions.
% Set to Inf to use all available acquisitions.
params.num_acquisitions_to_use = 100; 

% To speed up processing, you can downsample the data. A factor of 4-8
% is reasonable. This will reduce temporal resolution but speed up all steps.
% Set to 1 for no downsampling.
params.downsample_factor = 16; 

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
    fprintf('No background file selected. Proceeding without subtraction.\n');
end

% --- Automatically find matching profile files ---
fprintf('\n--- Automatically locating matching CSV profile files... ---\n');
[~, h5_name, ~] = fileparts(h5_file);
timestamp_match = regexp(h5_name, '\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', 'match');
if isempty(timestamp_match)
    error('Could not extract a valid timestamp (YYYY-MM-DD_HH-MM-SS) from the HDF5 filename.');
end
timestamp_str = timestamp_match{1};
fprintf('  Found timestamp in HDF5 filename: %s\n', timestamp_str);

freq_file = ['frequencies_' timestamp_str '.csv'];
delay_file = ['delays_' timestamp_str '.csv'];
phase_file = ['phases_' timestamp_str '.csv'];
freq_filepath = fullfile(h5_path, freq_file);
delay_filepath = fullfile(h5_path, delay_file);
phase_filepath = fullfile(h5_path, phase_file);

if ~exist(freq_filepath, 'file') || ~exist(delay_filepath, 'file') || ~exist(phase_filepath, 'file')
    error(['Could not automatically find all matching CSV files in the same directory as the HDF5 file.\n' ...
           'Please ensure these files are present:\n  - %s\n  - %s\n  - %s'], freq_file, delay_file, phase_file);
end
fprintf('Loaded profiles:\n  - Frequencies: %s\n  - Delays: %s\n  - Phases: %s\n', freq_file, delay_file, phase_file);

% --- Core Physical and Simulation Parameters ---
params.c = 343;
params.pmut_width_m = 0.020;
params.kerf_m = 0.001;
params.excitation_amplitude = 1e9;
params.filter_cutoff_hz = 70e3;
params.num_active_tx = 3;
params.grid_width_m = 0.150;
params.target_distance_m = 0.150;
params.grid_depth_range_m = 0.100;
params.grid_step_m = 0.002;
params.rho_admm = 7.0;
params.lambda_tv_reg = 1.5;
params.admm_tol = 1e-5;
params.admm_max_iter = 75;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;


%% ===== 2. LOAD & PRE-PROCESS EXPERIMENTAL DATA =====
fprintf('\n--- Step 1: Loading and Pre-processing Experimental Data ---\n');
tic;

% --- Load data and determine sampling rate ---
target_data_raw = h5read(data_filepath, '/echo_data_raw_adc');
exp_attrs = h5info(data_filepath, '/echo_data_raw_adc').Attributes;
timebase_attr = exp_attrs(strcmp({exp_attrs.Name}, 'timebase'));
timebase = double(timebase_attr.Value);

if isempty(timebase) || ~isnumeric(timebase) || timebase < 0
    error('Could not read a valid "timebase" attribute from HDF5 file.');
end

PICO_CLOCK_FREQ_HZ = 62.5e6;
if timebase < 3
    sample_interval_s = (2^timebase) / (PICO_CLOCK_FREQ_HZ * 2);
else
    sample_interval_s = (timebase - 2) / PICO_CLOCK_FREQ_HZ;
end
original_fs = 1 / sample_interval_s;
fprintf('  Detected original experimental sampling rate: %.2f MHz\n', original_fs / 1e6);

% --- Set final sampling rate based on downsampling factor ---
params.fs = original_fs / params.downsample_factor;
fprintf('  Effective sampling rate after downsampling by %dx: %.2f MHz\n', params.downsample_factor, params.fs / 1e6);

% --- Transpose data to correct orientation ---
target_data = target_data_raw';

% --- Perform Background Subtraction ---
if ~isempty(background_filepath)
    fprintf('  Performing background subtraction...\n');
    background_data_raw = h5read(background_filepath, '/echo_data_raw_adc');
    background_data = background_data_raw';
    if ~isequal(size(target_data), size(background_data))
        error('Target and Background data files have different dimensions after transpose.');
    end
    processed_data = double(target_data) - double(background_data);
else
    processed_data = double(target_data);
end

% --- Apply Low-Pass Filter (at original sampling rate to prevent aliasing) ---
fprintf('  Applying %.0f kHz low-pass filter...\n', params.filter_cutoff_hz / 1000);
try
    filtered_data = lowpass(processed_data, params.filter_cutoff_hz, original_fs);
catch ME
    if strcmp(ME.identifier, 'signal:lowpass:LicenseRequired')
        error('Signal Processing Toolbox license is required for the lowpass filter.');
    else
        rethrow(ME);
    end
end

% --- Downsample Data ---
if params.downsample_factor > 1
    fprintf('  Downsampling data by a factor of %d...\n', params.downsample_factor);
    filtered_data = downsample(filtered_data, params.downsample_factor);
end

% --- Load profiles and select subset if specified ---
fprintf('  Loading CSV profiles (skipping header row)...\n');
all_freqs = double(readmatrix(freq_filepath, 'NumHeaderLines', 1));
all_delays = double(readmatrix(delay_filepath, 'NumHeaderLines', 1));
all_phases = double(readmatrix(phase_filepath, 'NumHeaderLines', 1));

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

b_vector = final_data(:);

fprintf('Data loading and pre-processing complete. Time: %.2f seconds.\n', toc);


%% ===== 3. GENERATE H-MATRIX (DIGITAL TWIN) =====
fprintf('\n--- Step 2: Generating Digital Twin H-Matrix ---\n');
tic;
[H, imaging_grid] = generate_h_matrix_from_profiles(params);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

% --- Align H and b vectors ---
num_h_samples_per_acq = size(H, 1) / params.num_acquisitions;
num_b_samples_per_acq = size(final_data, 1);
if num_h_samples_per_acq ~= num_b_samples_per_acq
    fprintf('  Aligning H and b vectors (H samples: %d, b samples: %d)\n', ...
            num_h_samples_per_acq, num_b_samples_per_acq);
    min_samples = min(num_h_samples_per_acq, num_b_samples_per_acq);
    H_reshaped = reshape(full(H), [num_h_samples_per_acq, size(H, 2), params.num_acquisitions]);
    b_reshaped = reshape(b_vector, [num_b_samples_per_acq, 1, params.num_acquisitions]);
    H_aligned = H_reshaped(1:min_samples, :, :);
    b_aligned = b_reshaped(1:min_samples, :, :);
    H = sparse(reshape(H_aligned, [min_samples * params.num_acquisitions, size(H, 2)]));
    b_vector = b_aligned(:);
    fprintf('  Alignment complete. New vector length: %d\n', length(b_vector));
end


%% ===== 4. RECONSTRUCT IMAGE FROM REAL DATA =====
fprintf('\n--- Step 3: Reconstructing Image via ADMM ---\n');
tic;
reconstructed_image = run_admm_reconstruction(H, b_vector, zeros(size(imaging_grid.X_mesh)), params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);


%% ===== 5. ANALYZE AND SAVE RESULTS =====
fprintf('\n--- Step 4: Analyzing and Saving Results ---\n');
analyze_and_plot_reconstruction(reconstructed_image, imaging_grid, H, params, output_folder);
fprintf('\n\n=== RECONSTRUCTION COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix_from_profiles(config)
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);
    num_x_grid = 9; num_y_grid = 9;
    total_elements = num_x_grid * num_y_grid;
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
    for i = 1:num_active_tx
        distances = sum((element_centers - tx_desired_positions(i,:)).^2, 2);
        [~, tx_active_indices(i)] = min(distances);
    end
    tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    tx_enabled_matrix(tx_active_indices) = 1;
    TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
    [~, rx_active_index] = min(sum(element_centers.^2, 2));
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(rx_active_index) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Z_mesh', Z_mesh, 'points', grid_points);
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        tx_frequencies = repmat(config.exp_frequencies_hz(r_acq), 1, num_active_tx);
        delays_us = config.exp_delays_us(r_acq, :);
        phases_deg = config.exp_phases_deg(r_acq, :);
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k);
            t = 0:1/fs:duration;
            phase_rad = phases_deg(k) * pi / 180;
            signal_base = sin(2 * pi * tx_frequencies(k) * t + phase_rad);
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
        xdc_excitation(TxAperture, ones(1, num_active_tx));
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    disp('  Assembling H-matrix using interpolation...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), H = sparse(config.num_acquisitions, N_pixels); return; end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0, K_global_per_acq = 1; end
    total_rows = K_global_per_acq * config.num_acquisitions;
    H = spalloc(total_rows, N_pixels, round(sum(all_K_values) * N_pixels * 0.05));
    current_row_offset = 0;
    
    % *** NEW: Add waitbar for the slow assembly loop ***
    wb_asm = waitbar(0, 'Assembling H-matrix using interpolation...');
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned = interp1(t_current, all_h_data{r_acq}, t_common_axis, 'linear', 0);
            row_indices = current_row_offset + (1:K_global_per_acq);
            if max(row_indices) <= size(H, 1)
                H(row_indices, :) = h_aligned;
            end
        end
        current_row_offset = current_row_offset + K_global_per_acq;
        waitbar(r_acq / config.num_acquisitions, wb_asm); % Update waitbar
    end
    close(wb_asm); % Close waitbar
    
    xdc_free(TxAperture);
    xdc_free(RxAperture);
    field_end();
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
    fprintf('  Starting ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        x_prev = x_admm_img_iter;
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = admm_iteration(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            config.rho_admm, config.lambda_tv_reg, Hfun_pcg_admm, config);
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol
                fprintf('    ADMM converged at iteration %d.\n', k_admm);
                break;
            end
        end
        waitbar(k_admm / config.admm_max_iter, wb_admm, sprintf('ADMM Iteration %d/%d', k_admm, config.admm_max_iter));
    end
    close(wb_admm);
    reconstructed_image = x_admm_img_iter;
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [Dx, Dy] = difference_operators(imageSize)
   rows = imageSize(1); cols = imageSize(2);
   N_img_pixels = rows * cols;
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
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm = max(v_norm, eps);
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    u_new = u_old + opDx_tv(x_new) - z_new;
end

function analyze_and_plot_reconstruction(reconstructed_image, imaging_grid, H, params, output_folder)
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0;
    end
    fprintf('  Max Coherence of H-Matrix: %.4f\n', max_coherence);
    figure('Visible', 'on', 'Position', [100, 100, 700, 600]);
    imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Z_mesh(:,1)*1000, reconstructed_image);
    colormap(gray);
    colorbar;
    axis image;
    title(sprintf('Reconstructed Image from Real Data\n(Coherence: %.4f)', max_coherence));
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    saveas(gcf, fullfile(output_folder, 'final_reconstruction.png'));
    metrics = struct('max_coherence', max_coherence);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
    save(fullfile(output_folder, 'reconstruction_data.mat'), 'reconstructed_image', 'imaging_grid', 'params');
end
