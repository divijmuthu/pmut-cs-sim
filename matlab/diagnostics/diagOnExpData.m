% =========================================================================
% Combined Reconstruction & H-Matrix Diagnostic Script v1.0
%
% Description:
% This script provides a complete, end-to-end workflow. It first builds
% the H matrix using the definitive bistatic Field II model, and then
% immediately runs a comprehensive set of diagnostic analyses on that
% newly created matrix before proceeding to the final image reconstruction.
% =========================================================================

clear; clc; close all;

%% --- 1. User-Defined Parameters & File Selection ---

% --- Data Selection ---
ACQ_START_INDEX = 101;
ACQ_COUNT = 100;

% --- Physical & Simulation Constants ---
c = 343;
fs_picoscope = 41.67e6;

% --- SPEED & MEMORY OPTIMIZATION ---
downsample_factor = 40;
grid_step_m = 0.004;

fs = fs_picoscope / downsample_factor;

% --- PMUT Array Geometry ---
pmut_width_m = 0.020;
pmut_spacing_m = 0.020;
kerf_m = 0.0001;

% --- Imaging Grid Definition ---
grid_width_m = 0.150;
grid_depth_start_m = 0.250;
grid_depth_end_m = 0.350;

% --- Filter Parameters ---
filter_cutoff_hz = 80000;
filter_order = 4;

% --- Reconstruction Parameters ---
num_admm_iterations = 150;
rho_admm = 10;
lambda_tv_reg = 0.1;
max_iters_pcg = 100;
tol_pcg = 1e-6;
alpha_pcg_reg = 1e-4;
excitation_amplitude = 10000;

%% --- 2. Load and Pre-process Experimental Data ---
disp('--- Loading Experimental Data ---');

[h5_file, h5_path] = uigetfile('*.h5', 'Select the HDF5 Data File');
if isequal(h5_file, 0), disp('User selected Cancel'); return; end
h5_filepath = fullfile(h5_path, h5_file);
disp(['User selected H5 file: ', h5_filepath]);

[delays_file, delays_path] = uigetfile('*.csv', 'Select the Delays CSV File', h5_path);
if isequal(delays_file, 0), disp('User selected Cancel'); return; end
delays_filepath = fullfile(delays_path, delays_file);
disp(['User selected Delays file: ', delays_filepath]);

[freqs_file, freqs_path] = uigetfile('*.csv', 'Select the Frequencies CSV File', h5_path);
if isequal(freqs_file, 0), disp('User selected Cancel'); return; end
freqs_filepath = fullfile(freqs_path, freqs_file);
disp(['User selected Frequencies file: ', freqs_filepath]);

full_raw_adc_data = h5read(h5_filepath, '/echo_data_raw_adc');
full_delays_us = readmatrix(delays_filepath);
full_frequencies_hz = readmatrix(freqs_filepath);

if ACQ_START_INDEX + ACQ_COUNT - 1 > size(full_raw_adc_data, 2)
    error('Selected acquisition range exceeds the number of acquisitions in the file.');
end
acq_indices = ACQ_START_INDEX : (ACQ_START_INDEX + ACQ_COUNT - 1);
raw_adc_data = full_raw_adc_data(:, acq_indices);
delays_us = full_delays_us(acq_indices, :);
frequencies_hz = full_frequencies_hz(acq_indices);
num_acquisitions = size(raw_adc_data, 2);
fprintf('Selected %d acquisitions (from index %d to %d).\n', num_acquisitions, acq_indices(1), acq_indices(end));

max_val = max(abs(raw_adc_data), [], 'all');
b = double(raw_adc_data) / double(max_val);
b_resampled = resample(b, 1, downsample_factor);

disp('Applying low-pass filter to experimental data...');
nyquist = fs / 2;
Wn = filter_cutoff_hz / nyquist;
[b_filt, a_filt] = butter(filter_order, Wn, 'low');
b_filtered = filtfilt(b_filt, a_filt, b_resampled);

disp('Experimental data loaded and filtered.');

%% --- 3. Build H Matrix using Bistatic Simulation Recipe ---
disp('--- Building H Matrix with Field II ---');

field_init(-1);
set_field('fs', fs);
set_field('c', c);

tx_pos1 = [25e-3, 0, 0];
tx_pos2 = [-12.5e-3, 21.651e-3, 0];
tx_pos3 = [-12.5e-3, -21.651e-3, 0];
tx_desired_positions = [tx_pos1; tx_pos2; tx_pos3];
rx_pos = [0, 0, 0];
num_tx_pmuts = size(tx_desired_positions, 1);

num_x_grid = 9; num_y_grid = 9;
physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
center_offset_x = (num_x_grid-1)/2*(pmut_width_m+kerf_m); center_offset_y = (num_y_grid-1)/2*(pmut_width_m+kerf_m);
for iy = 1:num_y_grid
    y_pos_el = (iy-1)*(pmut_width_m+kerf_m) - center_offset_y;
    for ix = 1:num_x_grid
        x_pos_el = (ix-1)*(pmut_width_m+kerf_m) - center_offset_x;
        element_no = (iy-1)*num_x_grid + ix;
        physical_element_centers(element_no, :) = [x_pos_el, y_pos_el, 0];
    end
end

tx_active_indices = zeros(num_tx_pmuts, 1);
for i = 1:num_tx_pmuts
    distances = sum((physical_element_centers - tx_desired_positions(i, :)).^2, 2);
    [~, min_idx] = min(distances);
    tx_active_indices(i) = min_idx;
end
tx_active_indices = unique(tx_active_indices);
num_tx_active = length(tx_active_indices);

rx_distances = sum((physical_element_centers - rx_pos).^2, 2);
[~, rx_active_index] = min(rx_distances);

tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
tx_enabled_matrix(tx_active_indices) = 1;
tx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pmut_width_m, pmut_width_m, kerf_m, kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);

rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
rx_enabled_matrix(rx_active_index) = 1;
rx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pmut_width_m, pmut_width_m, kerf_m, kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);

f_start_chirp = 45e3;
f_end_chirp = 65e3;
burst_duration = 0.02e-3;
t_burst_vec = 0 : 1/fs : burst_duration;
synth_burst_base = chirp(t_burst_vec, f_start_chirp, t_burst_vec(end), f_end_chirp, 'linear');
synth_burst_windowed = synth_burst_base .* tukeywin(length(t_burst_vec), 0.25)';
impulse_response_waveform = synth_burst_windowed * excitation_amplitude;
xdc_impulse(tx_Aperture, impulse_response_waveform);
xdc_excitation(tx_Aperture, 1);
xdc_impulse(rx_Aperture, impulse_response_waveform);
xdc_excitation(rx_Aperture, 1);

x_coords_img = -grid_width_m/2 : grid_step_m : grid_width_m/2;
z_coords_img = grid_depth_start_m : grid_step_m : grid_depth_end_m;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
imageResolution = size(X_mesh);
grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
N_pixels = size(grid_points, 1);
fprintf('Imaging grid defined with %d pixels.\n', N_pixels);

all_h_data = cell(num_acquisitions, 1);
all_start_times = zeros(num_acquisitions, 1);
for r_acq = 1:num_acquisitions
    fprintf('  Simulating response for acquisition %d/%d...\n', r_acq, num_acquisitions);
    delays_for_sim = delays_us(r_acq, 1:num_tx_active) * 1e-6;
    xdc_focus_times(tx_Aperture, 0, delays_for_sim);
    [h_r, start_time_r] = calc_hhp(tx_Aperture, rx_Aperture, grid_points);
    all_h_data{r_acq} = h_r;
    all_start_times(r_acq) = start_time_r;
end

disp('Assembling full H matrix...');
K_per_acq = size(b_filtered, 1);
t_exp_axis = (0:K_per_acq-1) / fs;
H = zeros(K_per_acq * num_acquisitions, N_pixels);
for r_acq = 1:num_acquisitions
    h_r = all_h_data{r_acq};
    start_time_r = all_start_times(r_acq);
    t_sim_axis = start_time_r + (0:size(h_r, 1)-1) / fs;
    H_single_acq = zeros(K_per_acq, N_pixels);
    for px_col = 1:N_pixels
        H_single_acq(:, px_col) = interp1(t_sim_axis, h_r(:, px_col), t_exp_axis, 'linear', 0);
    end
    start_row = (r_acq-1)*K_per_acq + 1;
    end_row = r_acq*K_per_acq;
    H(start_row:end_row, :) = H_single_acq;
end
disp('H matrix assembled successfully.');

%% --- 4. H-MATRIX DIAGNOSTICS ---
disp('--- Running H-Matrix Diagnostics ---');
diag_output_folder = fullfile(h5_path, 'diagnostics_output');
if ~exist(diag_output_folder, 'dir'), mkdir(diag_output_folder); end
fprintf('Diagnostic figures will be saved to: %s\n', diag_output_folder);

% --- Diagnostic 4.1: System Energy Map ---
fprintf('  Calculating system energy map... ');
tic;
column_energies = vecnorm(H, 2, 1);
energy_map = reshape(column_energies, imageResolution);
toc;
figure('Name', 'Diagnostic 1: System Energy Map');
imagesc(x_coords_img*1000, z_coords_img*1000, energy_map);
title('System Energy Map (Sensitivity)');
xlabel('X Position (mm)'); ylabel('Z Position (mm)');
colorbar; colormap('hot'); axis image; set(gca, 'YDir', 'normal');
saveas(gcf, fullfile(diag_output_folder, 'diag_energy_map.png'));
fprintf('Saved energy map plot.\n');

% --- Diagnostic 4.2: Singular Value Decay ---
fprintf('  Calculating singular values... (This may take a moment)\n');
max_sv_cols = 500;
if N_pixels > max_sv_cols
    col_subset_indices = randperm(N_pixels, max_sv_cols);
    H_subset = H(:, col_subset_indices);
else
    H_subset = H;
end
tic;
S = svd(full(H_subset), 'econ');
toc;
figure('Name', 'Diagnostic 2: Singular Value Decay');
semilogy(S, 'LineWidth', 1.5);
title('Singular Value Decay of H Matrix (Subset)');
xlabel('Singular Value Index'); ylabel('Magnitude');
grid on; axis tight;
saveas(gcf, fullfile(diag_output_folder, 'diag_singular_values.png'));
fprintf('Saved singular value plot.\n\n');

%% --- 5. Automatic Time Alignment ---
disp('--- Automatically Aligning Data and Model ---');
b_vector_aligned = b_filtered(:);
center_pixel_idx = floor(N_pixels / 2);
h_sample = H(:, center_pixel_idx);
b_sample_for_xcorr = b_vector_aligned(1:K_per_acq);
[correlation, lags] = xcorr(b_sample_for_xcorr, h_sample(1:K_per_acq));
[~, max_lag_idx] = max(abs(correlation));
time_shift_samples = lags(max_lag_idx);
fprintf('  Detected time shift of %d samples (%.2f us).\n', time_shift_samples, time_shift_samples/fs*1e6);
b_vector_aligned = circshift(b_vector_aligned, -time_shift_samples);

%% --- 6. Reconstruction ---
fprintf('\n--- Starting Reconstruction ---\n');
A_matrix = H;
b_vector = b_vector_aligned;

% --- Fast Recon: PCG ---
fprintf('  Running PCG for a fast look...\n');
Afun = @(x, transp_flag) afun_for_pcg(x, A_matrix, transp_flag);
b_pcg = Afun(b_vector, 'transp');
pcg_system_fun = @(x) Afun(Afun(x,'notransp'),'transp') + alpha_pcg_reg * x;
tic;
[x_pcg_vec, flag_cg, ~, iter_cg] = pcg(pcg_system_fun, b_pcg, tol_pcg, max_iters_pcg);
runtime_pcg = toc;
fprintf('  PCG Finished: Flag=%d, Iterations=%d, Time=%.2f s\n', flag_cg, iter_cg, runtime_pcg);
x_pcg_img = reshape(real(x_pcg_vec), imageResolution);
figure;
imagesc(x_coords_img*100, z_coords_img*100, x_pcg_img);
title(sprintf('Fast Recon: PCG (%.2fs, %d iters)', runtime_pcg, iter_cg));
xlabel('X (cm)'); ylabel('Z (cm)'); colorbar; colormap('gray'); axis image; set(gca, 'YDir', 'normal');
drawnow;

% --- High-Quality Recon: ADMM-TV ---
fprintf('\n  Starting High-Quality ADMM-TV Reconstruction...\n');
H_norm_factor = max(abs(A_matrix(:)));
if H_norm_factor < eps, H_norm_factor = 1; end
A_admm = A_matrix ./ H_norm_factor;
b_admm_vec = b_vector(:) / H_norm_factor;
At_admm = transpose(A_admm);
Afun_admm = @(x) A_admm * x(:);
Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
[Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
x_admm_img_iter = zeros(imageResolution);
z_admm_grad_iter = zeros([prod(imageResolution) 2]);
u_admm_dual_iter = zeros([prod(imageResolution) 2]);
Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
figure;
tic;
for k_admm = 1:num_admm_iterations
    fprintf('  ADMM Iteration %d/%d\n', k_admm, num_admm_iterations);
    v_upd = z_admm_grad_iter - u_admm_dual_iter;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), 1e-8, 25, [], [], x_admm_img_iter(:));
    x_admm_img_iter = reshape(x_vec_new, imageResolution);
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm(v_norm < eps) = eps;
    shr = max(0, 1 - kap ./ v_norm);
    z_admm_grad_iter = v_z_upd .* shr;
    u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
    
    imagesc(x_coords_img*100, z_coords_img*100, real(x_admm_img_iter));
    title(sprintf('ADMM-TV Reconstruction (Iteration %d/%d)', k_admm, num_admm_iterations));
    xlabel('X (cm)'); ylabel('Z (cm)'); colorbar; colormap('gray'); axis image; set(gca, 'YDir', 'normal');
    drawnow;
end
toc;
disp('Reconstruction complete.');

%% --- Helper functions ---
function [Dx, Dy] = createDifferenceOperators(imageSize)
    rows = imageSize(1);
    cols = imageSize(2);
    N_img_pixels = rows * cols;
    Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
    last_col_indices_mask = false(N_img_pixels, 1);
    last_col_indices_mask((cols-1)*rows+1 : cols*rows) = true;
    Dx(last_col_indices_mask, :) = 0;
    Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
    last_row_indices_mask = false(N_img_pixels, 1);
    last_row_indices_mask(rows:rows:N_img_pixels) = true;
    Dy(last_row_indices_mask, :) = 0;
end

function y = afun_for_pcg(x, A, transp_flag)
   if strcmp(transp_flag,'transp')
      y = A'*x;
   else
      y = A*x;
   end
end
