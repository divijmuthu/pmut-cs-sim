% =========================================================================
% All-in-One Compressed Sensing Reconstruction v19 (Downsample & Filter)
%
% Description:
% This version fixes the filtering issue by first downsampling the raw
% data to a more reasonable rate, which allows for stable filter design.
% =========================================================================

clearvars;
clc;
close all;

% --- Configuration: Point to your data files ---
DELAY_PROFILE_CSV_FILE = 'delays_2025-07-07_17-03-00.csv';
MEASURED_DATA_H5_FILE = 'data_2025-07-07_17-03-00.h5';

% --- Time-Shift Correction Parameter ---
manual_time_shift_us = 375.12; 

% --- Initialize Field II ---
field_init(-1);

% --- Core Physical and Simulation Constants ---
c = 1540; fs = 125e6;
set_field('fs', fs); set_field('c', c);
fprintf('--- All-in-One Reconstruction Script v19 (Downsample & Filter) ---\n');

% --- Geometry and Simulation Parameters ---
pMUT_width_mm = 20; pMUT_spacing_mm = 25; kerf_mm = 0.1;
grid_width_mm = 150; grid_depth_start_mm = 250; grid_depth_end_mm = 350;
grid_step_mm = 8; 

excitation_amplitude = 500;
% --- Reconstruction Algorithm Parameters ---
maxItersCG_main = 200; tolCG_main = 1e-6; alpha_pcg_reg = 1e-4;
numItersADMM = 50; rho_admm = 10; lambda_tv_reg = 0.1;

% --- Convert mm to meters ---
pMUT_width = pMUT_width_mm/1000; pMUT_height = pMUT_width; kerf = kerf_mm/1000;
d_spacing = pMUT_spacing_mm/1000; grid_width = grid_width_mm/1000;
grid_depth_start = grid_depth_start_mm/1000; grid_depth_end = grid_depth_end_mm/1000;
grid_step = grid_step_mm/1000;

% --- Define Imaging Grid ---
x_coords_img = -grid_width/2 : grid_step : grid_width/2;
z_coords_img = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
Y_mesh = zeros(size(X_mesh));
N_pixels = numel(X_mesh);
imageResolution = [length(z_coords_img), length(x_coords_img)];
hydrophone_positions_img = [X_mesh(:), Y_mesh(:), Z_mesh(:)];
fprintf('Imaging Grid: %d pixels.\n', N_pixels);

% --- Define Separate Tx and Rx Apertures ---
fprintf('Defining separate Tx and Rx apertures...\n');
radius = d_spacing;
tx_pos1 = [radius, 0, 0]; tx_pos2 = [radius * cos(2*pi/3), radius * sin(2*pi/3), 0];
tx_pos3 = [radius * cos(4*pi/3), radius * sin(4*pi/3), 0];
tx_positions = [tx_pos1; tx_pos2; tx_pos3];
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

% --- Define Impulse Response ---
f_start_chirp = 10e3; f_end_chirp = 200e3; burst_duration = 0.02e-3;
t_burst_vec = 0 : 1/fs : burst_duration;
synth_burst_base = chirp(t_burst_vec, f_start_chirp, t_burst_vec(end), f_end_chirp, 'linear');
synth_burst_windowed = synth_burst_base .* tukeywin(length(t_burst_vec), 0.25)';
impulse_response_waveform = synth_burst_windowed * excitation_amplitude;
xdc_impulse(Tx_Aperture, impulse_response_waveform); xdc_excitation(Tx_Aperture, 1);

%% Load Experimental Data
fprintf('\n--- Loading Experimental Data ---\n');
if ~exist(DELAY_PROFILE_CSV_FILE, 'file'), error('Delay CSV not found: %s', DELAY_PROFILE_CSV_FILE); end
delay_profiles_us = readmatrix(DELAY_PROFILE_CSV_FILE);
if ~exist(MEASURED_DATA_H5_FILE, 'file'), error('HDF5 data file not found: %s', MEASURED_DATA_H5_FILE); end
echo_data_matrix = h5read(MEASURED_DATA_H5_FILE, '/echo_data')';
[R_acquisitions, total_samples_per_acq] = size(echo_data_matrix);
u_measured_signal_raw = double(reshape(echo_data_matrix', [], 1));
fprintf('Loaded data for %d acquisitions of %d samples each.\n', R_acquisitions, total_samples_per_acq);

%% --- NEW: Downsample and Apply Digital Low-Pass Filter ---
fprintf('\n--- Downsampling and Applying Digital Low-Pass Filter ---\n');
% Downsample data to a more manageable rate for filtering
downsample_factor = 50;
fs_new = fs / downsample_factor;
u_measured_downsampled = decimate(u_measured_signal_raw, downsample_factor);
fprintf('Downsampled signal from %.1f MHz to %.2f MHz.\n', fs/1e6, fs_new/1e6);

% Design a 4th-order Butterworth low-pass filter with a 60 kHz cutoff
cutoff_freq = 60e3; 
nyquist_freq = fs_new / 2;
filter_order = 4;
[b, a] = butter(filter_order, cutoff_freq / nyquist_freq, 'low');

% Apply the filter to the downsampled measured signal
u_measured_signal = filtfilt(b, a, u_measured_downsampled);
fprintf('Applied %d-order Butterworth filter with %.0f kHz cutoff.\n', filter_order, cutoff_freq/1e3);

% --- NEW: Diagnostic Plot for Filtering ---
figure(97); clf;
t_axis_us_raw = (0:total_samples_per_acq-1) / fs * 1e6;
t_axis_us_filt = (0:length(u_measured_signal(1:total_samples_per_acq/downsample_factor))-1) / fs_new * 1e6;
plot(t_axis_us_raw, echo_data_matrix(1, :), 'DisplayName', 'Raw Signal');
hold on;
plot(t_axis_us_filt, u_measured_signal(1:total_samples_per_acq/downsample_factor), 'r', 'LineWidth', 1.5, 'DisplayName', 'Filtered Signal');
hold off;
grid on; xlim([350 450]); % Zoom in on the echo region
legend; title('Filter Diagnostic'); xlabel('Time (\mus)');
fprintf('Displaying Figure 97 to show filter effect. Close to continue.\n');
uiwait(gcf);

%% Define Matrix-Free Operator Functions
% IMPORTANT: The operators must still work at the FULL sampling rate for the simulation
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
rms_measured = rms(u_measured_signal); % Use filtered signal for scaling
v_delta = zeros(N_pixels, 1); v_delta(floor(N_pixels/2)) = 1;
h_col_sample_raw = apply_H_operator(v_delta, all_hhp_data, all_start_times, R_acquisitions, total_samples_per_acq, fs, 0);
h_col_sample_filt = filtfilt(b, a, decimate(h_col_sample_raw, downsample_factor));
rms_simulated = rms(h_col_sample_filt);
scaling_factor = rms_measured / rms_simulated;
if isnan(scaling_factor) || isinf(scaling_factor) || scaling_factor == 0, scaling_factor = 1; end
fprintf('Scaling Factor: %.2e\n', scaling_factor);

Afun_final = @(v) filtfilt(b, a, decimate(scaling_factor * apply_H_operator(v, all_hhp_data, all_start_times, R_acquisitions, total_samples_per_acq, fs, manual_time_shift_us), downsample_factor));
Atfun_final = @(u) scaling_factor * apply_HT_operator(u, all_hhp_data, all_start_times, R_acquisitions, total_samples_per_acq, fs, N_pixels, manual_time_shift_us, b, a, downsample_factor);

%% --- Reconstruction Algorithms ---
fprintf('\n--- Starting Reconstruction ---\n');
b_vector = u_measured_signal;

% --- 1. Least Norm (PCG) with Regularization ---
fprintf('  Calculating Least Norm (PCG)...\n'); tic;
AAtfun_pcg_reg = @(y_vec) Afun_final(Atfun_final(y_vec)) + alpha_pcg_reg * y_vec;
[y_sol, ~, ~, ~] = pcg(AAtfun_pcg_reg, b_vector(:), tolCG_main, maxItersCG_main);
x_pcg_img_vec = Atfun_final(y_sol);
runtime_pcg = toc; fprintf('  PCG finished in %.2f s.\n', runtime_pcg);
x_pcg2D = reshape(real(x_pcg_img_vec), imageResolution);

% --- 2. ADMM with TV regularization ---
fprintf('  Calculating ADMM-TV Solution...\n'); tic;
[Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
x_admm_img_iter = zeros(imageResolution);
z_admm_grad_iter = zeros([prod(imageResolution) 2]);
u_admm_dual_iter = zeros([prod(imageResolution) 2]);
Atb_vector = Atfun_final(b_vector);
for k_admm = 1:numItersADMM
    fprintf('    ADMM Iteration %d/%d\n', k_admm, numItersADMM);
    v_upd = z_admm_grad_iter - u_admm_dual_iter;
    opDtx_v = reshape(Dx_sparse' * v_upd(:, 1) + Dy_sparse' * v_upd(:, 2), imageResolution);
    bb_upd = Atb_vector + rho_admm * opDtx_v(:);
    AtAfun_op = @(x) Afun_final(Afun_final(x(:)));
    opDtDx_op = @(x) reshape(Dx_sparse'*(Dx_sparse*x(:)) + Dy_sparse'*(Dy_sparse*x(:)), size(x));
    Hfun_pcg_admm = @(x) AtAfun_op(x) + rho_admm * opDtDx_op(reshape(x,imageResolution));
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
set(gcf, 'Position', [100, 100, 1000, 500], 'Color', 'w');
normalize = @(img) (img - min(img(:))) / (max(img(:)) - min(img(:)));

subplot(1, 2, 1);
imagesc(x_coords_img, z_coords_img, normalize(x_pcg2D));
axis image; colormap gray; colorbar; set(gca, 'YDir', 'normal');
title(sprintf('Least Norm (PCG + Regularization)\nRuntime: %.2fs', runtime_pcg));
xlabel('Lateral (m)'); ylabel('Axial (m)');

subplot(1, 2, 2);
imagesc(x_coords_img, z_coords_img, normalize(x_admm2D));
axis image; colormap gray; colorbar; set(gca, 'YDir', 'normal');
title(sprintf('ADMM-TV (Iters=%d)\nRuntime: %.2fs', numItersADMM, runtime_ADMM));
xlabel('Lateral (m)');

sgtitle('Comparison of Reconstructions (Matrix-Free, Regularized)');
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
        t_sim_axis = start_time_r + (0:(K_current-1))/fs + time_shift_s; % Apply shift
        t_target_axis = (0:(K_acq-1))/fs;
        echo_aligned = interp1(t_sim_axis, echo_r_sim, t_target_axis, 'linear', 0);
        start_idx = (r_acq-1) * K_acq + 1;
        end_idx = r_acq * K_acq;
        y(start_idx:end_idx) = echo_aligned;
    end
end

function v = apply_HT_operator(u, all_hhp_data, all_start_times, R, K_acq, fs, N_pixels, time_shift_us, b, a, ds_factor)
    v = zeros(N_pixels, 1);
    time_shift_s = time_shift_us / 1e6;
    fs_new = fs / ds_factor;
    for r_acq = 1:R
        hhp_r = all_hhp_data{r_acq};
        if isempty(hhp_r), continue; end
        
        start_idx_ds = (r_acq-1) * (K_acq/ds_factor) + 1;
        end_idx_ds = r_acq * (K_acq/ds_factor);
        u_r = u(start_idx_ds:end_idx_ds);
        
        start_time_r = all_start_times(r_acq);
        K_current = size(hhp_r, 1);
        t_sim_axis = start_time_r + (0:(K_current-1))/fs + time_shift_s;
        t_target_axis = (0:(K_acq-1))/fs;
        hhp_aligned_r = interp1(t_sim_axis, hhp_r, t_target_axis, 'linear', 0);
        
        hhp_filt = filtfilt(b, a, decimate(hhp_aligned_r, ds_factor));
        
        v = v + hhp_filt' * u_r;
    end
end

function [Dx, Dy] = createDifferenceOperators(imageSize)
    rows = imageSize(1); cols = imageSize(2); N = rows * cols;
    Dx = spdiags([-ones(N, 1), ones(N, 1)], [0, rows], N, N);
    Dx(N-rows+1:N, :) = 0;
    Dy = spdiags([-ones(N, 1), ones(N, 1)], [0, 1], N, N);
    Dy(rows:rows:N, :) = 0;
end
