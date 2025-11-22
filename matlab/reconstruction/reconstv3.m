% =========================================================================
% All-in-One Compressed Sensing Reconstruction v29 (Definitive)
%
% Description:
% This is the definitive reconstruction script, incorporating all learned
% corrections: 3Tx/1Rx pitch-catch model, aggressive downsampling to make
% the computation feasible, and the critical time-shift and amplitude
% scaling to align the simulation with the physical data.
% =========================================================================

clearvars;
clc;
close all;

% --- Configuration: Point to your data files ---
DELAY_PROFILE_CSV_FILE = 'delays_2025-07-07_17-03-00.csv';
MEASURED_DATA_H5_FILE = 'data_2025-07-07_17-03-00.h5';

% --- Time-Shift Correction Parameter ---
manual_time_shift_us = 375.12; % From previous diagnostic

% --- Core Physical and Simulation Constants ---
c = 1540;                   % Speed of Sound [m/s]
fs_physical = 125e6;        % The TRUE sampling rate of the collected data
downsample_factor = 62;
fs = fs_physical / downsample_factor; % Effective sampling rate is now ~2 MHz

% --- Initialize Field II ---
field_init(-1);
set_field('fs', fs);
set_field('c', c);
fprintf('--- Reconstruction with Downsampling ---\n');
fprintf('Original Fs: %.1f MHz | Downsample Factor: %d | New Fs: %.2f MHz\n\n', fs_physical/1e6, downsample_factor, fs/1e6);

% --- Geometry and Simulation Parameters ---
pMUT_width_mm = 20; pMUT_spacing_mm = 25; kerf_mm = 0.1;
grid_width_mm = 150; 
grid_depth_start_mm = 100;  % Start imaging at 10cm
grid_depth_end_mm = 200;    % End imaging at 20cm
grid_step_mm = 2;           % Using a finer grid step for better resolution

excitation_amplitude = 500;
% --- Reconstruction Algorithm Parameters ---
maxItersCG_main = 1000; tolCG_main = 1e-8;
numItersADMM = 100; rho_admm = 10; lambda_tv_reg = 0.1;

% --- Convert mm to meters ---
pMUT_width = pMUT_width_mm/1000; pMUT_height = pMUT_width; kerf = kerf_mm/1000;
d_spacing = pMUT_spacing_mm/1000; grid_width = grid_width_mm/1000;
grid_depth_start = grid_depth_start_mm/1000; grid_depth_end = grid_depth_end_mm/1000;
grid_step = grid_step_mm/1000;

% --- Define Imaging Grid ---
x_coords_img = -grid_width/2 : grid_step : grid_width/2;
z_coords_img = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
N_pixels = numel(X_mesh);
imageResolution = [length(z_coords_img), length(x_coords_img)];
hydrophone_positions_img = [X_mesh(:), zeros(N_pixels, 1), Z_mesh(:)];
fprintf('Imaging Grid: %d pixels, Depth: %.0fmm - %.0fmm\n', N_pixels, grid_depth_start_mm, grid_depth_end_mm);

% --- Define Separate Tx and Rx Apertures using Grid Mapping ---
fprintf('Defining separate Tx and Rx apertures...\n');
radius = d_spacing;
tx_positions = [radius, 0, 0; radius*cos(2*pi/3), radius*sin(2*pi/3), 0; radius*cos(4*pi/3), radius*sin(4*pi/3), 0];
rx_position = [0, 0, 0];
num_x_grid = 15; num_y_grid = 15;
physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
element_no_grid_map = 0;
center_offset_x = (num_x_grid - 1) / 2 * (pMUT_width + kerf);
center_offset_y = (num_y_grid - 1) / 2 * (pMUT_height + kerf);
for iy = 1:num_y_grid, y_pos_el = (iy - 1) * (pMUT_height + kerf) - center_offset_y;
for ix = 1:num_x_grid, x_pos_el = (ix - 1) * (pMUT_width + kerf) - center_offset_x;
element_no_grid_map = element_no_grid_map + 1;
physical_element_centers(element_no_grid_map, :) = [x_pos_el, y_pos_el, 0]; end, end
tx_indices_linear = zeros(size(tx_positions,1), 1);
for i = 1:size(tx_positions,1), [~, tx_indices_linear(i)] = min(sum((physical_element_centers - tx_positions(i,:)).^2, 2)); end
tx_indices_linear = unique(tx_indices_linear); num_active_tx = length(tx_indices_linear);
tx_enabled_matrix = zeros(num_y_grid, num_x_grid); tx_enabled_matrix(tx_indices_linear) = 1;
Tx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, tx_enabled_matrix, 1, 1, [0 0 0]);
[~, rx_index_linear] = min(sum((physical_element_centers - rx_position).^2, 2));
rx_enabled_matrix = zeros(num_y_grid, num_x_grid); rx_enabled_matrix(rx_index_linear) = 1;
Rx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, rx_enabled_matrix, 1, 1, [0 0 0]);
num_active = num_active_tx;
fprintf('Successfully mapped %d Tx and 1 Rx element.\n', num_active);

% --- Define Impulse Response ---
f_start_chirp = 10e3; f_end_chirp = 200e3; burst_duration = 0.02e-3;
t_burst_vec = 0 : 1/fs : burst_duration;
synth_burst_base = chirp(t_burst_vec, f_start_chirp, t_burst_vec(end), f_end_chirp, 'linear');
synth_burst_windowed = synth_burst_base .* tukeywin(length(t_burst_vec), 0.25)';
impulse_response_waveform = synth_burst_windowed * excitation_amplitude;
xdc_impulse(Tx_Aperture, impulse_response_waveform); xdc_excitation(Tx_Aperture, 1);

%% Load and Process Experimental Data
fprintf('\n--- Loading and Processing Experimental Data ---\n');
if ~exist(DELAY_PROFILE_CSV_FILE, 'file'), error('Delay CSV not found!'); end
delay_profiles_us = readmatrix(DELAY_PROFILE_CSV_FILE);
if ~exist(MEASURED_DATA_H5_FILE, 'file'), error('HDF5 data file not found!'); end
info = h5info(MEASURED_DATA_H5_FILE);
dataset_name = info.Datasets(1).Name;
echo_data_matrix_raw = h5read(MEASURED_DATA_H5_FILE, ['/' dataset_name])';
[R_acquisitions, ~] = size(echo_data_matrix_raw);

% --- Downsample the Measured Data ---
temp_decimated = decimate(double(echo_data_matrix_raw(1,:)), downsample_factor);
u_measured_downsampled = zeros(R_acquisitions, length(temp_decimated));
u_measured_downsampled(1, :) = temp_decimated;
for i = 2:R_acquisitions
    u_measured_downsampled(i, :) = decimate(double(echo_data_matrix_raw(i,:)), downsample_factor);
end
u_measured_signal = reshape(u_measured_downsampled', [], 1);
total_samples_per_acq_new = size(u_measured_downsampled, 2);
fprintf('Data downsampled. New samples per acquisition: %d\n', total_samples_per_acq_new);

%% Define Matrix-Free Operator Functions
fprintf('\n--- Defining Matrix-Free Operators ---\n');
all_hhp_data = cell(R_acquisitions, 1);
all_start_times = zeros(R_acquisitions, 1);
for r_acq = 1:R_acquisitions
    fprintf('Pre-calculating impulse response for Acq %d/%d...\n', r_acq, R_acquisitions);
    delay_vector_s = delay_profiles_us(r_acq, :) / 1e6;
    xdc_focus_times(Tx_Aperture, 0, delay_vector_s);
    [hhp_r, start_time_r] = calc_hhp(Tx_Aperture, Rx_Aperture, hydrophone_positions_img);
    all_hhp_data{r_acq} = hhp_r;
    all_start_times(r_acq) = start_time_r;
end

%% Automatic Scaling
fprintf('\n--- Automatically Scaling Simulation ---\n');
rms_measured = rms(u_measured_signal);
v_delta = zeros(N_pixels, 1); v_delta(floor(N_pixels/2)) = 1;
h_col_sample = apply_H_operator(v_delta, all_hhp_data, all_start_times, R_acquisitions, total_samples_per_acq_new, fs, 0);
rms_simulated = rms(h_col_sample);
if abs(rms_simulated) > 1e-20
    scaling_factor = rms_measured / rms_simulated;
else
    scaling_factor = 1;
end
fprintf('Applied scaling factor of %.2e to H matrix.\n', scaling_factor);

Afun_final = @(v) scaling_factor * apply_H_operator(v, all_hhp_data, all_start_times, R_acquisitions, total_samples_per_acq_new, fs, manual_time_shift_us);
Atfun_final = @(u) scaling_factor * apply_HT_operator(u, all_hhp_data, all_start_times, R_acquisitions, total_samples_per_acq_new, fs, N_pixels, manual_time_shift_us);

%% --- Reconstruction Algorithms ---
fprintf('\n--- Starting Reconstruction ---\n');
b_vector = u_measured_signal;

% --- 1. Least Norm (PCG) ---
fprintf('  Calculating Least Norm (PCG)...\n'); tic;
AAtfun_pcg = @(y_vec) Afun_final(Atfun_final(y_vec));
[y_sol, ~, ~, ~] = pcg(AAtfun_pcg, b_vector(:), tolCG_main, maxItersCG_main);
x_pcg_img_vec = Atfun_final(y_sol);
runtime_pcg = toc; fprintf('  PCG finished in %.2f s.\n', runtime_pcg);
x_pcg2D = reshape(real(x_pcg_img_vec), imageResolution);

% --- 2. Least Norm (Pseudo-inverse) ---
fprintf('  Calculating Least Norm (Pinv)...\n'); tic;
% Note: Pinv is not feasible for matrix-free, so we build a smaller H
H_for_pinv = zeros(total_samples_per_acq_new * R_acquisitions, N_pixels);
for i = 1:N_pixels
    v_temp = zeros(N_pixels, 1); v_temp(i) = 1;
    H_for_pinv(:, i) = Afun_final(v_temp);
end
try
    x_pinv_img_vec = pinv(H_for_pinv) * b_vector;
    runtime_pinv = toc; fprintf('  Pinv finished in %.2f s.\n', runtime_pinv);
    x_pinv2D = reshape(real(x_pinv_img_vec), imageResolution);
catch ME
    warning('Pinv failed, likely due to matrix size. Skipping. Error: %s');
    x_pinv2D = zeros(imageResolution); runtime_pinv = toc;
end

% --- 3. ADMM with TV regularization ---
fprintf('  Calculating ADMM-TV Solution...\n'); tic;
[Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
x_admm_img_iter = zeros(imageResolution);
z_admm_grad_iter = zeros([prod(imageResolution) 2]);
u_admm_dual_iter = zeros([prod(imageResolution) 2]);
Atb_vector = Atfun_final(b_vector);
AtAfun_op = @(x_img) reshape(Atfun_final(Afun_final(x_img(:))), imageResolution);
opDtDx_op = @(x_img) reshape(Dx_sparse'*(Dx_sparse*x_img(:)) + Dy_sparse'*(Dy_sparse*x_img(:)), size(x_img));
Hfun_pcg_admm = @(x_vec) reshape(AtAfun_op(reshape(x_vec, imageResolution)) + rho_admm * opDtDx_op(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
for k_admm = 1:numItersADMM
    fprintf('    ADMM Iteration %d/%d\n', k_admm, numItersADMM);
    v_upd = z_admm_grad_iter - u_admm_dual_iter;
    opDtx_v = reshape(Dx_sparse' * v_upd(:, 1) + Dy_sparse' * v_upd(:, 2), imageResolution);
    bb_upd = Atb_vector + rho_admm * opDtx_v(:);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), 1e-5, 25, [], [], x_admm_img_iter(:));
    x_admm_img_iter = reshape(x_vec_new, imageResolution);
    kap = lambda_tv_reg / rho_admm;
    opDx_v = [Dx_sparse * x_admm_img_iter(:), Dy_sparse * x_admm_img_iter(:)];
    v_z_upd = opDx_v + u_admm_dual_iter;
    v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm(v_norm < eps) = eps;
    shr = max(0, 1 - kap ./ v_norm);
    z_admm_grad_iter = v_z_upd .* shr;
    u_admm_dual_iter = u_admm_dual_iter + opDx_v - z_admm_grad_iter;
end
runtime_ADMM = toc; fprintf('  ADMM-TV finished in %.2f s.\n', runtime_ADMM);
x_admm2D = reshape(real(x_admm_img_iter), imageResolution);

%% --- Plot Final Reconstructions ---
figure(99); clf;
set(gcf, 'Position', [100, 100, 1400, 500], 'Color', 'w');
normalize = @(img) (img - min(img(:))) / (max(img(:)) - min(img(:)));

subplot(1, 3, 1);
imagesc(x_coords_img*1000, z_coords_img*1000, normalize(x_pcg2D));
axis image; colormap gray; colorbar; set(gca, 'YDir', 'normal');
title(sprintf('Least Norm (PCG)\nRuntime: %.2fs', runtime_pcg));
xlabel('Lateral (mm)'); ylabel('Axial (mm)');

subplot(1, 3, 2);
imagesc(x_coords_img*1000, z_coords_img*1000, normalize(x_pinv2D));
axis image; colormap gray; colorbar; set(gca, 'YDir', 'normal');
title(sprintf('Least Norm (Pinv)\nRuntime: %.2fs', runtime_pinv));
xlabel('Lateral (mm)');

subplot(1, 3, 3);
imagesc(x_coords_img*1000, z_coords_img*1000, normalize(x_admm2D));
axis image; colormap gray; colorbar; set(gca, 'YDir', 'normal');
title(sprintf('ADMM-TV (Iters=%d)\nRuntime: %.2fs', numItersADMM, runtime_ADMM));
xlabel('Lateral (m)');

sgtitle('Comparison of Reconstructions (Downsampled to ~2 MHz)');
fprintf('\n--- Reconstruction Complete ---\n');

%% End Field II Simulation
if exist('field_end', 'file') == 2, field_end; disp('Field II ended.'); end

%% --- Matrix-Free Operator Functions ---
function y = apply_H_operator(v, all_hhp_data, all_start_times, R, K_acq, fs, time_shift_us)
    y = zeros(R * K_acq, 1);
    time_shift_s = time_shift_us / 1e6;
    for r_acq = 1:R
        hhp_r = all_hhp_data{r_acq};
        if isempty(hhp_r), continue; end
        echo_r_sim = hhp_r * v;
        start_time_r = all_start_times(r_acq);
        K_current = size(hhp_r, 1);
        t_sim_axis = start_time_r + (0:(K_current-1))/fs + time_shift_s;
        t_target_axis = (0:(K_acq-1))/fs;
        echo_aligned = interp1(t_sim_axis, echo_r_sim, t_target_axis, 'linear', 0);
        start_idx = (r_acq-1) * K_acq + 1;
        end_idx = r_acq * K_acq;
        y(start_idx:end_idx) = echo_aligned;
    end
end

function v = apply_HT_operator(u, all_hhp_data, all_start_times, R, K_acq, fs, N_pixels, time_shift_us)
    v = zeros(N_pixels, 1);
    time_shift_s = time_shift_us / 1e6;
    for r_acq = 1:R
        hhp_r = all_hhp_data{r_acq};
        if isempty(hhp_r), continue; end
        
        start_idx = (r_acq-1) * K_acq + 1;
        end_idx = r_acq * K_acq;
        u_r = u(start_idx:end_idx);
        
        start_time_r = all_start_times(r_acq);
        K_current = size(hhp_r, 1);
        t_sim_axis = start_time_r + (0:(K_current-1))/fs + time_shift_s;
        t_target_axis = (0:(K_acq-1))/fs;
        hhp_aligned_r = interp1(t_sim_axis, hhp_r, t_target_axis, 'linear', 0);
        
        v = v + hhp_aligned_r' * u_r;
    end
end

function [Dx, Dy] = createDifferenceOperators(imageSize)
    rows = imageSize(1); cols = imageSize(2); N = rows * cols;
    Dx = spdiags([-ones(N, 1), ones(N, 1)], [0, rows], N, N);
    Dx(N-rows+1:N, :) = 0;
    Dy = spdiags([-ones(N, 1), ones(N, 1)], [0, 1], N, N);
    Dy(rows:rows:N, :) = 0;
end
