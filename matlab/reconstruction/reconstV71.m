% =========================================================================
% FINAL REAL DATA RECONSTRUCTION SCRIPT (v4.8 - Final Tuning)
%
% Description:
% This version is designed for the final, crucial step of hyperparameter
% tuning. The key ADMM parameters (`lambda_tv_reg` and `rho_admm`) have
% been moved to the top for easy experimentation. Adjusting these is the
% standard method for improving image quality when working with noisy,
% real-world data.
%
% v4.8 Improvements:
% - Prominently featured ADMM parameters at the top for easy tuning.
% - Added detailed comments explaining what each parameter does.
% =========================================================================

clear; clc; close all;

%% ===== 1. CONFIGURATION AND FILE SELECTION =====
fprintf('=== pMUT FINAL RECONSTRUCTION (v4.8) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- ADMM HYPERPARAMETER TUNING ---
% This is the most important section for improving final image quality.
%
% lambda_tv_reg: The "Clean-Up" or "Denoising" parameter.
%   - Higher values (e.g., 10, 25, 50) make the algorithm more aggressive
%     at removing noise and artifacts, resulting in a smoother image.
%   - Lower values (e.g., 1.2) trust the data more, which can result
%     in a sharper but noisier image.
%   - RECOMMENDATION: Start by increasing this value significantly.
params.lambda_tv_reg = 25.0; % (Original v1.5 value was 1.2)

% rho_admm: The "Step Size" or "Penalty" parameter.
%   - This affects the convergence speed and stability. It's best to start
%     by tuning lambda_tv_reg first.
params.rho_admm = 6.73; % (Original v1.5 value)


% --- Performance & Debugging Settings ---
params.num_acquisitions_to_use = 120; 
params.downsample_factor = 4;
params.admm_max_iter = 75;

% --- Select Experimental Data Files ---
[h5_file, h5_path] = uigetfile('*.h5', 'Select the TARGET HDF5 data file');
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
params.pmut_width_m = 0.002; 
params.kerf_m = 0.018;       
params.excitation_amplitude = 1e9;
params.filter_cutoff_hz = 70e3;
params.num_active_tx = 3;

% --- Imaging Grid Parameters (Top-Down XY Plane) ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.003;

% --- Time Gating Parameters ---
params.gate_start_time_s = (params.target_height_m / params.c) * 0.8;
params.gate_duration_s = 3.0e-3;

% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
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


%% ===== 3. GENERATE H-MATRIX COMPONENTS =====
fprintf('\n--- Step 2: Generating H-Matrix Components for Forward Operator ---\n');
tic;
[h_components, imaging_grid, h_norm_factor] = generate_h_matrix_components(params);
fprintf('H-component generation complete. Time: %.2f seconds.\n', toc);
fprintf('  Estimated H-matrix normalization factor: %.2e\n', h_norm_factor);

% --- Align and Prune Data ---
num_h_samples_per_acq = size(h_components{1}, 1);
num_b_samples_per_acq = size(final_data, 1);
if num_h_samples_per_acq ~= num_b_samples_per_acq
    fprintf('  Aligning H-components and b vector...\n');
    min_samples = min(num_h_samples_per_acq, num_b_samples_per_acq);
    for i = 1:params.num_acquisitions
        h_components{i} = h_components{i}(1:min_samples, :);
    end
    b_reshaped = reshape(b_vector, [num_b_samples_per_acq, params.num_acquisitions]);
    b_aligned = b_reshaped(1:min_samples, :);
    b_vector = b_aligned(:);
    num_samples_per_acq = min_samples;
else
    num_samples_per_acq = num_b_samples_per_acq;
end

fprintf('  Applying time gate and pruning zero-rows...\n');
time_axis = (0:num_samples_per_acq-1)' / params.fs;
gate_mask = (time_axis >= params.gate_start_time_s) & ...
            (time_axis < (params.gate_start_time_s + params.gate_duration_s));

for i = 1:params.num_acquisitions
    h_components{i} = h_components{i}(gate_mask, :);
end
b_reshaped = reshape(b_vector, [num_samples_per_acq, params.num_acquisitions]);
b_pruned = b_reshaped(gate_mask, :);
b_final = b_pruned(:);

fprintf('  Pruning complete. Kept %d samples per acquisition.\n', size(b_pruned, 1));


%% ===== 4. RECONSTRUCT IMAGE FROM REAL DATA =====
fprintf('\n--- Step 3: Reconstructing Image via ADMM ---\n');
tic;
reconstructed_image = run_admm_reconstruction(h_components, b_final, zeros(size(imaging_grid.X_mesh)), h_norm_factor, params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);


%% ===== 5. ANALYZE AND SAVE RESULTS =====
fprintf('\n--- Step 4: Analyzing and Saving Results ---\n');
analyze_and_plot_reconstruction(reconstructed_image, imaging_grid, output_folder);
fprintf('\n\n=== RECONSTRUCTION COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

function [all_h_aligned, imaging_grid, h_norm_factor] = generate_h_matrix_components(config)
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
    s = 0.047;
    tx_desired_positions = [ 0, (sqrt(3)/3)*s, 0; -s/2, -(sqrt(3)/6)*s, 0; s/2, -(sqrt(3)/6)*s, 0 ];
    tx_active_indices = zeros(1, num_active_tx);
    for i = 1:num_active_tx, [~, tx_active_indices(i)] = min(sum((element_centers - tx_desired_positions(i,:)).^2, 2)); end
    
    FullTxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, ones(num_y_grid,num_x_grid), 1, 1, [0 0 0]);
    
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
        if mod(r_acq, 10) == 0, fprintf('  Generating H-matrix acquisition %d/%d (%.0f%%)...\n', r_acq, config.num_acquisitions, (r_acq/config.num_acquisitions)*100); end
        
        center_freq = config.exp_frequencies_hz(r_acq);
        duration = 3 / center_freq;
        t = 0:1/fs:duration;
        impulse_response = sin(2*pi*center_freq*t) .* tukeywin(length(t), 0.25)';
        xdc_impulse(FullTxAperture, impulse_response * config.excitation_amplitude);

        apodization_weights = zeros(1, total_elements);
        apodization_weights(tx_active_indices) = config.exp_cycles(r_acq, :) / max(config.exp_cycles(:));
        xdc_apodization(FullTxAperture, 0, apodization_weights);
        
        delays_and_phases_s = zeros(1, total_elements);
        delays_s = config.exp_delays_us(r_acq, :) * 1e-6;
        phases_s = (config.exp_phases_deg(r_acq, :) / 360) / center_freq;
        delays_and_phases_s(tx_active_indices) = delays_s + phases_s;
        xdc_focus_times(FullTxAperture, 0, delays_and_phases_s);

        [h_r, start_time_r] = calc_hhp(FullTxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r; all_start_times(r_acq) = start_time_r; all_K_values(r_acq) = size(h_r, 1);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    disp('  Aligning H-matrix components to common time axis...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), error('No valid H-matrix components were generated.'); end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices)); max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    
    all_h_aligned = cell(config.num_acquisitions, 1);
    max_h_val = 0;
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned_acq = interp1(t_current, all_h_data{r_acq}, t_common_axis, 'linear', 0);
            all_h_aligned{r_acq} = h_aligned_acq;
            max_h_val = max(max_h_val, max(abs(h_aligned_acq(:))));
        else
            all_h_aligned{r_acq} = zeros(length(t_common_axis), N_pixels);
        end
    end
    h_norm_factor = max_h_val;
    if h_norm_factor < eps, h_norm_factor = 1; end
    
    xdc_free(FullTxAperture); xdc_free(RxAperture); field_end();
end

function reconstructed_image = run_admm_reconstruction(h_components, b_vector, scene_matrix_size, h_norm_factor, config)
    imageResolution = size(scene_matrix_size);
    
    [Afun, Atfun, AtAfun] = operator_setup(h_components, imageResolution, h_norm_factor, config);
    
    b_admm_vec = b_vector(:) / h_norm_factor;
    
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
    
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    fprintf('  Starting ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 5) == 0 || k_admm == config.admm_max_iter, fprintf('  ADMM iteration %d/%d (%.0f%%)...\n', k_admm, config.admm_max_iter, (k_admm/config.admm_max_iter)*100); end
        x_prev = x_admm_img_iter;
        
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
        
        [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter);
        x_admm_img_iter = reshape(x_vec_new, imageResolution);
        
        kap = config.lambda_tv_reg / config.rho_admm;
        v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
        v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm = max(v_norm, eps);
        shr = max(0, 1 - kap ./ v_norm);
        z_admm_grad_iter = v_z_upd .* shr;
        u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
        
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol, fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        end
        waitbar(k_admm / config.admm_max_iter, wb_admm, sprintf('ADMM Iteration %d/%d', k_admm, config.admm_max_iter));
    end
    close(wb_admm);
    reconstructed_image = x_admm_img_iter;
end

function [Afun, Atfun, AtAfun] = operator_setup(h_components, imageResolution, h_norm_factor, config)
    num_acqs = config.num_acquisitions;
    num_pixels = prod(imageResolution);
    num_samples_per_acq = size(h_components{1}, 1);

    function y = forward_op(x_img)
        x_vec = x_img(:);
        y = zeros(num_samples_per_acq * num_acqs, 1);
        for i = 1:num_acqs
            start_idx = (i-1) * num_samples_per_acq + 1;
            end_idx = i * num_samples_per_acq;
            y(start_idx:end_idx) = (h_components{i} * x_vec) / h_norm_factor;
        end
    end

    function x_img = transpose_op(y_vec)
        y_reshaped = reshape(y_vec, [num_samples_per_acq, num_acqs]);
        x_vec = zeros(num_pixels, 1);
        for i = 1:num_acqs
            x_vec = x_vec + (h_components{i}' * y_reshaped(:, i)) / h_norm_factor;
        end
        x_img = reshape(x_vec, imageResolution);
    end

    Afun = @(x) forward_op(x);
    Atfun = @(y) transpose_op(y);
    AtAfun = @(x) transpose_op(forward_op(x));
end

function [Dx, Dy] = difference_operators(imageSize)
   rows = imageSize(1); cols = imageSize(2); N_img_pixels = rows * cols;
   Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
   Dx( (cols-1)*rows+1 : cols*rows , :) = 0;
   Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
   Dy( rows:rows:N_img_pixels , :) = 0;
end

function analyze_and_plot_reconstruction(reconstructed_image, imaging_grid, output_folder)
    figure('Visible', 'on', 'Position', [100, 100, 700, 600]);
    imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Y_mesh(:,1)*1000, reconstructed_image);
    colormap(gray); colorbar; axis image;
    title('Reconstructed Image from Real Data');
    xlabel('X Position (mm)'); ylabel('Y Position (mm)');
    saveas(gcf, fullfile(output_folder, 'final_reconstruction.png'));
    save(fullfile(output_folder, 'reconstruction_data.mat'), 'reconstructed_image', 'imaging_grid');
end
