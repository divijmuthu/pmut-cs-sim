% Advanced_pMUT_CS_Recon_SynthChirp_v2.m
%
% Description:
% Version 2 of the advanced compressed sensing (CS) ultrasound simulation.
%
% --- KEY FIXES & IMPROVEMENTS ---
% 1.  Syntax Correction: All non-breaking space characters that can cause
%     "Invalid expression" errors have been removed.
% 2.  Broadband Transducer Model: The transducer's impulse response is now
%     set to a delta function ('impulseResponse = [1];'). This models a very
%     broadband pMUT that can faithfully transmit the full 10-200 kHz chirp
%     excitation. This is intended to fix the "all gray image" issue by
%     ensuring a strong signal is transmitted, leading to a robust H matrix.
%
% --- EXCITATION SOURCE ---
% This version GENERATES a SYNTHETIC chirp waveform (10-200 kHz) internally.
%
% Date: 2025-05-29
%

clearvars; clc; close all;

% --- Initialize Field II ---
field_init(-1);

% --- Simulation Parameters ---
fc = 5.8e4;       % Nominal Center Frequency [Hz] (Used for lambda reference)
c = 1540;         % Speed of Sound [m/s]
lambda = c/fc;    % Wavelength [m]
fs = 2e6;         % Simulation Sampling Frequency [Hz]. Set to a high fixed value for chirp.

% Set Field Parameters
set_field('fs', fs);
set_field('c', c);

% --- Coded Acquisition Parameters ---
R_acquisitions = 20; % Number of unique coded acquisitions (Increase for better quality)

% --- Define Impulse Response & Global Excitation ---
% IMPROVEMENT: Use a delta function for the impulse response to model a
% transducer that is extremely broadband. This ensures the chirp from
% ele_waveform is not filtered out by a narrow impulse response, which was
% a likely cause of weak signals and gray reconstructions.
impulseResponse = [1];
fprintf('Set impulse response to a delta function for a broadband transducer model.\n');
excitationPulse = 1;

% --- pMUT Element Parameters & Array Definition ---
pMUT_width = 1.5*lambda;
pMUT_height = pMUT_width;
kerf = 0.5*lambda;
d_spacing = 5*lambda;
pos1 = [d_spacing,0,0];
pos2 = [d_spacing*cos(2*pi/3), d_spacing*sin(2*pi/3),0];
pos3 = [d_spacing*cos(4*pi/3), d_spacing*sin(4*pi/3),0];
desired_positions = [pos1; pos2; pos3];
num_active_intended = size(desired_positions, 1);
array_span_x = 2.5*d_spacing;
array_span_y = 2.5*d_spacing;
element_width = pMUT_width;
element_height = pMUT_height;
num_x = ceil(array_span_x / (element_width + kerf));
num_y = ceil(array_span_y / (element_height + kerf));
if mod(num_x, 2) == 0; num_x = num_x + 1; end
if mod(num_y, 2) == 0; num_y = num_y + 1; end
num_elements_total = num_x * num_y;
physical_element_centers = zeros(num_elements_total, 3);
element_no_grid = 0;
center_offset_x = (num_x - 1)/2 * (element_width + kerf);
center_offset_y = (num_y - 1)/2 * (element_height + kerf);
for iy = 1:num_y
    y_pos_el = (iy-1)*(element_height + kerf) - center_offset_y;
    for ix = 1:num_x
        x_pos_el = (ix-1)*(element_width + kerf) - center_offset_x;
        element_no_grid = element_no_grid + 1;
        physical_element_centers(element_no_grid, :) = [x_pos_el, y_pos_el, 0];
    end
end
active_indices_linear = zeros(num_active_intended, 1);
for i = 1:num_active_intended
    distances = sqrt(sum((physical_element_centers - desired_positions(i,:)).^2, 2));
    [~, min_idx] = min(distances);
    active_indices_linear(i) = min_idx;
end
active_indices_linear = unique(active_indices_linear);
num_active = length(active_indices_linear);
if num_active ~= num_active_intended
    warning('Could not map to %d unique elements (found %d).', num_active_intended, num_active);
    if num_active == 0; error('No active elements found!'); end
end
enabled = zeros(num_y, num_x);
[row_indices, col_indices] = ind2sub([num_y, num_x], active_indices_linear);
for i = 1:num_active
    enabled(row_indices(i), col_indices(i)) = 1;
end
focus_point = [0 0 100];
sub_x = 1;
sub_y = 1;
Th = xdc_2d_array(num_x, num_y, element_width, element_height, kerf, kerf, enabled, sub_x, sub_y, focus_point);
xdc_impulse(Th, impulseResponse);
xdc_excitation(Th, excitationPulse);

% --- Imaging Grid Definition ---
grid_width = 20*lambda;
grid_depth_start = 10*lambda;
grid_depth_end = 40*lambda;
grid_step = lambda/2;
x_grid = -grid_width/2 : grid_step : grid_width/2;
z_grid = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_grid, z_grid);
Y_mesh = zeros(size(X_mesh));
N_pixels = numel(X_mesh);
fprintf('Imaging Grid: %d pixels.\n', N_pixels);

%% Generate Measurement Matrix H from Multiple Coded Acquisitions (SYNTHETIC CHIRP VERSION)
fprintf('\n--- Generating Measurement Matrix H from %d Coded Acquisitions ---\n', R_acquisitions);
fprintf('--- Using SYNTHETIC Chirp (10-200 kHz) ---\n');

hydrophone_positions = [X_mesh(:), Y_mesh(:), Z_mesh(:)];

% =================== Generate Synthetic Chirp ===================
f_start = 10e3;  % 10 kHz
f_end = 200e3; % 200 kHz
chirp_duration = 1.5e-3; % 1.5 ms duration (Adjust as needed)
t_chirp = 0 : 1/fs : chirp_duration;
synth_chirp = chirp(t_chirp, f_start, t_chirp(end), f_end, 'linear');
synth_chirp_windowed = synth_chirp .* tukeywin(length(t_chirp), 0.25)';
base_chirp_waveform_sim = synth_chirp_windowed(:)'; % Ensure it's a row vector

% Optional: Plot the synthetic chirp to verify
figure(99); clf;
plot(t_chirp * 1e3, base_chirp_waveform_sim);
title('Generated Synthetic Chirp Pulse'); xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
saveas(gcf, 'Advanced_pMUT_Synthetic_Chirp.png');
fprintf('Generated synthetic chirp pulse (10-200 kHz, %g ms duration).\n', chirp_duration*1e3);
% ===============================================================

all_hhp_data = cell(R_acquisitions, 1);
all_start_times = zeros(R_acquisitions, 1);
all_K_values = zeros(R_acquisitions, 1);
rng('default');

for r = 1:R_acquisitions
    fprintf('Calculating Acquisition %d/%d...\n', r, R_acquisitions);
    for p_idx = 1:num_active
        current_rand_phase = 2 * pi * rand();
        analytic_signal = hilbert(base_chirp_waveform_sim);
        phase_shifted_analytic = analytic_signal * exp(1i * current_rand_phase);
        current_pmut_waveform = real(phase_shifted_analytic);
        ele_waveform(Th, p_idx, current_pmut_waveform);
    end
    [hhp_r, start_time_r] = calc_hhp(Th, Th, hydrophone_positions);
    all_hhp_data{r} = hhp_r;
    all_start_times(r) = start_time_r;
    all_K_values(r) = size(hhp_r, 1);
    if ~isempty(hhp_r)
        fprintf('    Max abs value in hhp_r: %g\n', max(abs(hhp_r(:))));
    else
        fprintf('    WARNING: hhp_r is empty for acquisition %d!\n', r);
    end
end

%% Assemble Full H Matrix with Interpolation and Alignment
fprintf('\n--- Assembling Full H Matrix ---\n');
all_end_times = all_start_times + (all_K_values - 1)/fs;
all_end_times(all_K_values == 0) = all_start_times(all_K_values == 0);
min_global_start_time = min(all_start_times);
max_global_end_time = max(all_end_times);
t_common = min_global_start_time : 1/fs : max_global_end_time;
K_global = length(t_common);
fprintf('Global Time Window: K_global_per_acq = %d samples.\n', K_global);
H_final = zeros(K_global * R_acquisitions, N_pixels);
current_row_offset = 0;
for r = 1:R_acquisitions
    hhp_current = all_hhp_data{r};
    start_time_current = all_start_times(r);
    K_current = all_K_values(r);
    if K_current == 0
        current_row_offset = current_row_offset + K_global;
        continue;
    end
    t_current_acq = start_time_current + (0:(K_current-1))/fs;
    hhp_aligned_r = zeros(K_global, N_pixels);
    for px_col = 1:N_pixels
        hhp_aligned_r(:, px_col) = interp1(t_current_acq, hhp_current(:,px_col), t_common, 'linear', 0);
    end
    H_final(current_row_offset + (1:K_global), :) = hhp_aligned_r;
    current_row_offset = current_row_offset + K_global;
end
H = H_final;
M = size(H, 1);
N = N_pixels;
fprintf('Final Measurement Matrix H: %d rows (M) x %d columns (N).\n', M, N);

%% Generate Scene (v) - Point targets
fprintf('\n--- Generating Scene (v) ---\n');
scene_xmin = min(x_grid); scene_xmax = max(x_grid);
scene_zmin = min(z_grid); scene_zmax = max(z_grid); zrange = scene_zmax - scene_zmin;
targets.x = [scene_xmin*0.5, 0, scene_xmax*0.5];
targets.y = zeros(size(targets.x));
targets.z = [scene_zmin + 0.2*zrange, scene_zmin + 0.5*zrange, scene_zmin + 0.8*zrange];
amplitudes = [1, 0.9, 1.1]';
scene = zeros(length(z_grid), length(x_grid));
for i = 1:length(targets.x)
    [~, ix] = min(abs(x_grid - targets.x(i)));
    [~, iz] = min(abs(z_grid - targets.z(i)));
    if ix > 0 && ix <= length(x_grid) && iz > 0 && iz <= length(z_grid)
        scene(iz, ix) = amplitudes(i);
    end
end
v = scene(:);

%% Image formation model: u = H*v + n
fprintf('\n--- Simulating Measurement Acquisition (u = Hv + n) ---\n');
Hv = H * v;
electronic_SNR_db = 25;
electronic_SNR = 10^(electronic_SNR_db/10);
signal_power = mean(Hv(:).^2);
if isnan(signal_power) || signal_power < 1e-20
    noise_sigma = 1e-7;
    actual_SNR_db = -Inf;
else
    noise_variance = signal_power / electronic_SNR;
    noise_sigma = sqrt(noise_variance);
    actual_SNR_db = 10*log10(signal_power / noise_variance);
end
n_noise = noise_sigma * randn(size(Hv));
u = Hv + n_noise;
fprintf('Added Gaussian noise. Target SNR: %.1f dB. Actual SNR: %.1f dB.\n', electronic_SNR_db, actual_SNR_db);

%% Reconstruction Algorithms
fprintf('\n--- Starting Reconstruction ---\n');
A = H;
At = transpose(A);
b_meas = u;
I_true = scene;
imageResolution = size(I_true);
Afun = @(x_vec) A*x_vec;
Atfun = @(y_vec) At*y_vec;
norm_I_true_vec = I_true(:);
if max(abs(norm_I_true_vec)) > 0
    norm_I_true_vec = norm_I_true_vec ./ max(abs(norm_I_true_vec));
end

% --- PCG Reconstruction ---
tic;
fprintf('  Calculating Least Norm (PCG) Solution...\n');
AAtfun = @(x_intermediate_cg_vec) Afun(Atfun(x_intermediate_cg_vec));
maxItersCG = 200;
tolCG = noise_sigma;
if M==0 || N==0 || isempty(b_meas) || norm(b_meas)==0
    x_pcg_img_vec = zeros(N,1);
    iter_cg = 0;
else
    [x_intermediate_cg_vec, flag_cg] = pcg(AAtfun, b_meas(:), tolCG, maxItersCG);
    if flag_cg ~= 0; warning('PCG did not converge (Flag=%d).', flag_cg); end
    x_pcg_img_vec = Atfun(x_intermediate_cg_vec);
end
x_pcg_img_vec = real(x_pcg_img_vec);
min_val_pcg = min(x_pcg_img_vec(:));
max_val_pcg = max(x_pcg_img_vec(:));
if (max_val_pcg - min_val_pcg) > eps
    x_pcg_norm_vec = (x_pcg_img_vec - min_val_pcg) / (max_val_pcg - min_val_pcg);
else
    x_pcg_norm_vec = zeros(N,1);
end
MSE_pcg = mean((x_pcg_norm_vec - norm_I_true_vec).^2);
PSNR_pcg = 10*log10(1/MSE_pcg);
runtime_pcg = toc;
fprintf('  PCG: PSNR = %g dB, Time = %g s\n', PSNR_pcg, runtime_pcg);

% --- Pinv Reconstruction ---
tic;
fprintf('  Calculating Least Norm (Pinv) Solution...\n');
tolPinv = 1e-6;
x_pinv_img_vec = zeros(N, 1);
pinv_method_desc = 'Skipped';
if M>0 && N>0 && ~isempty(A) && ~isempty(b_meas)
    try
        if M < N
            AAt = A*At;
            AAt_pinv = pinv(AAt, tolPinv);
            x_pinv_img_vec = At * (AAt_pinv * b_meas(:));
            pinv_method_desc = 'At*pinv(A*At)*b';
        else
            A_pinv = pinv(A, tolPinv);
            x_pinv_img_vec = A_pinv * b_meas(:);
            pinv_method_desc = 'pinv(A)*b';
        end
    catch ME
        fprintf('Pinv failed:%s\n', ME.message);
        pinv_method_desc = 'Failed';
    end
end
x_pinv_img_vec = real(x_pinv_img_vec);
min_val_pinv = min(x_pinv_img_vec(:));
max_val_pinv = max(x_pinv_img_vec(:));
if (max_val_pinv - min_val_pinv) > eps
    x_pinv_norm_vec = (x_pinv_img_vec - min_val_pinv) / (max_val_pinv - min_val_pinv);
else
    x_pinv_norm_vec = zeros(N,1);
end
MSE_pinv = mean((x_pinv_norm_vec - norm_I_true_vec).^2);
PSNR_pinv = 10*log10(1/MSE_pinv);
runtime_pinv = toc;
fprintf('  Pinv (%s): PSNR = %g dB, Time = %g s\n', pinv_method_desc, PSNR_pinv, runtime_pinv);


% --- ADMM Reconstruction ---
fprintf('  Calculating ADMM-TV Solution...\n');
H_norm_factor = max(abs(H(:)));
if H_norm_factor < eps; H_norm_factor = 1; end
A_admm = H ./ H_norm_factor;
At_admm = transpose(A_admm);
b_admm_vec = b_meas(:) / H_norm_factor;
noise_sigma_admm = noise_sigma / H_norm_factor;
Afun_admm = @(x_vec) A_admm*x_vec;
Atfun_admm_img = @(y_vec) reshape(At_admm*y_vec, imageResolution);
AtAfun_admm_img = @(x_img) Atfun_admm_img(Afun_admm(x_img(:)));
[Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
opDx = @(x_img) [Dx_sparse*x_img(:), Dy_sparse*x_img(:)];
opDtx = @(v_vecs) reshape(Dx_sparse'*v_vecs(:,1) + Dy_sparse'*v_vecs(:,2), imageResolution);
opDtDx = @(x_img) opDtx(opDx(x_img));
numItersADMM = 30;
rho = 1.0;
lambda_tv = 0.01;
x_admm_img = zeros(imageResolution);
z_admm_grad_vecs = zeros([prod(imageResolution), 2]);
u_admm_dual_vecs = zeros([prod(imageResolution), 2]);
Hfun_admm_pcg = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec,imageResolution)) + rho.*opDtDx(reshape(x_vec,imageResolution)), [prod(imageResolution), 1]);
PSNR_admm_iters = zeros([numItersADMM, 1]);
tic;
for k_admm = 1:numItersADMM
    rhs_admm_pcg_img = Atfun_admm_img(b_admm_vec) + rho*opDtx(z_admm_grad_vecs - u_admm_dual_vecs);
    maxItersCG_admm = 25;
    tolCG_admm = 1e-5;
    if norm(rhs_admm_pcg_img(:)) == 0
        x_admm_img_vec_new = zeros(prod(imageResolution), 1);
    else
        [x_admm_img_vec_new] = pcg(Hfun_admm_pcg, rhs_admm_pcg_img(:), tolCG_admm, maxItersCG_admm, [], [], x_admm_img(:));
    end
    x_admm_img = reshape(x_admm_img_vec_new, imageResolution);
    kappa_tv = lambda_tv/rho;
    v_grad_update = opDx(x_admm_img) + u_admm_dual_vecs;
    v_grad_norm = sqrt(sum(v_grad_update.^2, 2));
    v_grad_norm(v_grad_norm < eps) = eps;
    shrinkage_factor = max(0, 1 - kappa_tv ./ v_grad_norm);
    z_admm_grad_vecs = v_grad_update .* shrinkage_factor;
    u_admm_dual_vecs = u_admm_dual_vecs + opDx(x_admm_img) - z_admm_grad_vecs;
    x_admm_scaled_vec = real(x_admm_img(:));
    min_x = min(x_admm_scaled_vec);
    max_x = max(x_admm_scaled_vec);
    if (max_x-min_x) > eps
        x_admm_scaled_vec = (x_admm_scaled_vec - min_x) / (max_x - min_x);
    else
        x_admm_scaled_vec = zeros(size(x_admm_scaled_vec));
    end
    MSE_admm = mean((x_admm_scaled_vec - norm_I_true_vec).^2);
    PSNR_admm_iters(k_admm) = 10*log10(1/MSE_admm);
end
runtime_ADMM_total = toc;
fprintf('  ADMM-TV complete: Final PSNR = %g dB, Time = %gs\n', PSNR_admm_iters(end), runtime_ADMM_total);

% Code to plot final results (not included in the compact version for brevity)

%% End Field II Simulation
if exist('field_end','file') == 2
    field_end;
    disp('Field II ended.');
end

%% Helper function for TV operators
function [Dx, Dy] = createDifferenceOperators(imageSize)
    rows = imageSize(1);
    cols = imageSize(2);
    N = rows*cols;
    % Difference in y-direction (vertical)
    Dy = spdiags([-ones(N,1), ones(N,1)], [0, 1], N, N);
    % Boundary condition for last row
    Dy(rows:rows:N, :) = 0;
    % Difference in x-direction (horizontal)
    Dx = spdiags([-ones(N,1), ones(N,1)], [0, rows], N, N);
    % Boundary condition for last column
    Dx(N-rows+1:N, :) = 0;
end

