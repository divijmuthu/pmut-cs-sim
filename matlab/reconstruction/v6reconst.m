% =========================================================================
% Ultrasound Reconstruction Script v2.5 (Final Fixes)
%
% Description:
% This version provides the final fixes to produce a valid image.
%
% Key Improvements:
% 1. Fixes the data-type crash in the ADMM solver permanently.
% 2. Removes all pre-trigger data to solve the time-alignment issue and
%    prevent the reconstruction from locking onto the EMI artifact.
% =========================================================================
clearvars;
clc;
close all;

%% --- 1. Configuration ---
fprintf('--- 1. Configuring Reconstruction ---\n');

% --- Input Files ---
DELAY_PROFILE_CSV_FILE = 'delays_2025-07-11_14-38-54.csv';
MEASURED_DATA_H5_FILE  = 'data_2025-07-11_14-38-54.h5';

% --- Reconstruction Parameters (SPEED & MEMORY) ---
DOWNSAMPLING_FACTOR = 40;
DATA_TRUNCATE_MS = 3;
IMAGE_GRID_SIZE = 64;

% --- Reconstruction Parameters (Algorithm) ---
LAMBDA_TV = 0.5;
ADMM_ITERATIONS = 50;
ADMM_RHO = 1.0;
PCG_ITERATIONS = 20;

% --- Physical & Simulation Parameters ---
IMAGING_DEPTH_M = 0.25;
IMAGING_WIDTH_M = 0.25;

%% --- 2. Load and Pre-process Experimental Data ---
fprintf('\n--- 2. Loading and Pre-processing Data ---\n');

% --- Load Delays ---
if ~exist(DELAY_PROFILE_CSV_FILE, 'file'), error('Delay CSV not found!'); end
delay_profiles_us = readmatrix(DELAY_PROFILE_CSV_FILE);
[R_acquisitions, ~] = size(delay_profiles_us);

% --- Load H5 Data ---
if ~exist(MEASURED_DATA_H5_FILE, 'file'), error('HDF5 data file not found!'); end
info = h5info(MEASURED_DATA_H5_FILE);
dataset_name = info.Datasets(1).Name;
echo_data_matrix_raw = h5read(MEASURED_DATA_H5_FILE, ['/' dataset_name])';

% --- ADC-to-Voltage conversion block ---
try
    adc_max = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'max_adc');
    adc_min = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'min_adc');
    v_max = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'max_volts');
    v_min = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'min_volts');
    m = (v_max - v_min) / (adc_max - adc_min);
    c_offset = v_max - m * adc_max;
    echo_data_matrix_volts = double(echo_data_matrix_raw) * m + c_offset;
catch
    warning('Could not read voltage attributes. Using raw data as voltage.');
    echo_data_matrix_volts = double(echo_data_matrix_raw);
end
clear echo_data_matrix_raw;

% --- Calculate True Sampling Rate ---
[~, total_samples_raw] = size(echo_data_matrix_volts);
total_duration_s = 10e-3;
fs_raw = total_samples_raw / total_duration_s;
fprintf('Original sampling rate: %.2f MHz\n', fs_raw/1e6);

% --- Lowpass Filter ---
fc_cutoff = 80e3;
[b_filter, a_filter] = butter(4, fc_cutoff / (fs_raw/2), 'low');
echo_data_matrix_filtered = zeros(size(echo_data_matrix_volts));
for i = 1:R_acquisitions
    echo_data_matrix_filtered(i, :) = filtfilt(b_filter, a_filter, echo_data_matrix_volts(i, :));
end
clear echo_data_matrix_volts;

% *** FIX: Remove pre-trigger data to solve alignment issues ***
trigger_sample_raw = round(0.1 * total_samples_raw);
echo_data_post_trigger = echo_data_matrix_filtered(:, trigger_sample_raw:end);
clear echo_data_matrix_filtered;

% --- Truncate data ---
samples_to_keep_raw = round(DATA_TRUNCATE_MS / 1000 * fs_raw);
if size(echo_data_post_trigger, 2) < samples_to_keep_raw
    samples_to_keep_raw = size(echo_data_post_trigger, 2);
    fprintf('Warning: Truncation duration is longer than available post-trigger data.\n');
end
echo_data_truncated = echo_data_post_trigger(:, 1:samples_to_keep_raw);
fprintf('Data truncated to first %d ms post-trigger.\n', DATA_TRUNCATE_MS);
clear echo_data_post_trigger;

% --- Simple Decimation for Downsampling ---
fs_new = fs_raw / DOWNSAMPLING_FACTOR;
echo_data_downsampled = echo_data_truncated(:, 1:DOWNSAMPLING_FACTOR:end);
[~, total_samples_new] = size(echo_data_downsampled);
b = echo_data_downsampled(:);
fprintf('Downsampled to %.2f MHz. New data size: %d acquisitions x %d samples.\n', fs_new/1e6, R_acquisitions, total_samples_new);
clear echo_data_downsampled echo_data_truncated;

%% --- 3. Setup Field II Simulation Environment ---
fprintf('\n--- 3. Setting up Field II Simulation ---\n');
field_init(-1);
set_field('fs', fs_new);
set_field('c', 1540);
% (Transducer definition code is unchanged)
pMUT_width = 20/1000; pMUT_height = pMUT_width; kerf = 0.1/1000;
d_spacing = 25 / 1000; radius = d_spacing;
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
tx_indices_linear = unique(tx_indices_linear);
tx_enabled_matrix = zeros(num_y_grid, num_x_grid); tx_enabled_matrix(tx_indices_linear) = 1;
Tx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, tx_enabled_matrix, 1, 1, [0 0 0]);
[~, rx_index_linear] = min(sum((physical_element_centers - rx_position).^2, 2));
rx_enabled_matrix = zeros(num_y_grid, num_x_grid); rx_enabled_matrix(rx_index_linear) = 1;
Rx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, rx_enabled_matrix, 1, 1, [0 0 0]);
f0 = 55e3;
impulse_response = sin(2*pi*f0*(0:1/fs_new:2/f0));
xdc_impulse(Tx_Aperture, impulse_response);
xdc_excitation(Tx_Aperture, 1);

%% --- 4. Build System Matrix A ---
fprintf('\n--- 4. Building System Matrix A (This will be much faster now) ---\n');
N = IMAGE_GRID_SIZE;
x_range = linspace(-IMAGING_WIDTH_M/2, IMAGING_WIDTH_M/2, N);
z_range = linspace(0, IMAGING_DEPTH_M, N);
[X, Z] = meshgrid(x_range, z_range);
image_pixels = [X(:), zeros(N*N, 1), Z(:)];
num_pixels = N*N;

% --- Simplified time axis, starting from t=0 (the trigger) ---
time_axis_new = (0:total_samples_new-1) / fs_new;

A = zeros(length(b), num_pixels, 'single');
h_wait = waitbar(0, 'Building System Matrix A...');
for i = 1:R_acquisitions
    waitbar(i/R_acquisitions, h_wait, sprintf('Processing Acquisition %d/%d', i, R_acquisitions));
    delay_vector_s = delay_profiles_us(i, :) / 1e6;
    xdc_focus_times(Tx_Aperture, 0, delay_vector_s);
    [hhp_sim, start_time_sim] = calc_hhp(Tx_Aperture, Rx_Aperture, image_pixels);
    for p = 1:num_pixels
        if isscalar(start_time_sim), current_start_time = start_time_sim; else, current_start_time = start_time_sim(p); end
        t_sim_axis = current_start_time + (0:(size(hhp_sim,1)-1))/fs_new;
        u_sim_raw = interp1(t_sim_axis, hhp_sim(:,p), time_axis_new, 'linear', 0);
        u_sim_raw(~isfinite(u_sim_raw)) = 0;
        u_sim_filtered = filtfilt(b_filter, a_filter, u_sim_raw);
        start_idx = (i-1)*total_samples_new + 1;
        end_idx = i*total_samples_new;
        A(start_idx:end_idx, p) = single(u_sim_filtered(:));
    end
end
close(h_wait);
fprintf('System Matrix A built successfully. Size: %d x %d\n', size(A,1), size(A,2));

%% --- 5. Solve using Multiple Algorithms ---
fprintf('\n--- 5. Solving with Multiple Algorithms ---\n');

% --- Pre-calculate for speed ---
A_T = A';
A_TA = A_T * A;
A_Tb = A_T * b;
clear A; A_T; % Free memory

% --- Algorithm 1: Pseudoinverse (Fastest, Baseline) ---
fprintf('Solving with Pseudoinverse...\n');
x_pinv = pinv(A_TA) * A_Tb;
fprintf('Done.\n');

% --- Display Pseudoinverse result immediately ---
figure('Name', 'Immediate Result: Pseudoinverse');
img_pinv_temp = reshape(x_pinv, N, N);
imagesc(x_range*100, z_range*100, img_pinv_temp);
colormap('gray'); colorbar; axis image;
title('Immediate Preview: Pseudoinverse');
xlabel('Lateral (cm)'); ylabel('Axial (cm)');
drawnow;

% --- Algorithm 2: Preconditioned Conjugate Gradient (PCG) ---
fprintf('Solving with PCG...\n');
[x_pcg, ~] = pcg(A_TA, A_Tb, 1e-6, PCG_ITERATIONS);
fprintf('Done.\n');

% --- Algorithm 3: ADMM with TV Regularization ---
fprintf('Solving with ADMM...\n');
[Grad, Div] = get_div_grad_operators(N, N);
x_admm = zeros(num_pixels, 1, 'single');
z = zeros(size(Grad(x_admm)), 'single');
u = zeros(size(z), 'single');
h_wait_admm = waitbar(0, 'Running ADMM Solver...');
for k = 1:ADMM_ITERATIONS
    waitbar(k/ADMM_ITERATIONS, h_wait_admm, sprintf('Iteration %d/%d', k, ADMM_ITERATIONS));
    x_admm = (A_TA + ADMM_RHO * (Grad' * Grad)) \ (A_Tb + ADMM_RHO * Grad' * (z - u));
    Grad_x = Grad(x_admm);
    z = soft_threshold(Grad_x + u, LAMBDA_TV / ADMM_RHO);
    u = u + Grad_x - z;
end
close(h_wait_admm);
fprintf('Done.\n');

%% --- 6. Display Final Comparison ---
fprintf('\n--- 6. Displaying Final Comparison ---\n');
img_pinv = reshape(x_pinv, N, N);
img_pcg = reshape(x_pcg, N, N);
img_admm = reshape(x_admm, N, N);

% --- Standardize color scale for better comparison ---
min_val = min([img_pinv(:); img_pcg(:); img_admm(:)]);
max_val = max([img_pinv(:); img_pcg(:); img_admm(:)]);
clims = [min_val max_val];

figure('Position', [100, 100, 1500, 500], 'Name', 'Final Algorithm Comparison');

subplot(1,3,1);
imagesc(x_range*100, z_range*100, img_pinv, clims);
colormap('gray'); colorbar; axis image;
xlabel('Lateral (cm)'); ylabel('Axial (cm)');
title('1. Pseudoinverse (Least Squares)');

subplot(1,3,2);
imagesc(x_range*100, z_range*100, img_pcg, clims);
colormap('gray'); colorbar; axis image;
xlabel('Lateral (cm)');
title(sprintf('2. PCG (%d Iterations)', PCG_ITERATIONS));

subplot(1,3,3);
imagesc(x_range*100, z_range*100, img_admm, clims);
colormap('gray'); colorbar; axis image;
xlabel('Lateral (cm)');
title(sprintf('3. ADMM-TV (%d Iterations)', ADMM_ITERATIONS));

sgtitle('Reconstruction Algorithm Comparison');

field_end;
disp('Script finished.');

%% --- Helper Functions ---
function [Grad, Div] = get_div_grad_operators(n, m)
    % Create sparse operators as double, then cast the input vector
    D_x = spdiags([-ones(n,1), ones(n,1)], [0, 1], n, n);
    D_y = spdiags([-ones(m,1), ones(m,1)], [0, 1], m, m);
    D_x(n,n) = 0; D_y(m,m) = 0;
    % *** FIX: Cast the input vector 'I' to double for the multiplication ***
    Grad = @(I) [D_y * double(I); double(I) * D_x'];
    Div = @(P) -( (D_y' * P(1:n*m,:)) + (P(n*m+1:end,:)*D_x)' );
end

function z = soft_threshold(x, kappa)
    z = max(0, 1 - kappa./sqrt(sum(x.^2, 2))) .* x;
end
