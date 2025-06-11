% Advanced_pMUT_CS_Recon_SynthChirp_v5.m
%
% Description:
% Version 3 of the advanced compressed sensing (CS) ultrasound simulation.
% This version aims to:
% 1.  Improve PCG speed by reducing chirp duration and default R_acquisitions.
% 2.  More closely adhere to the structure and plotting of the original
%     Ensemble.m script for familiarity.
% 3.  Retain the core improvements: synthetic chirp excitation, broadband
%     transducer model, and robust H-matrix assembly for the 3-pMUT system.
%
% --- EXCITATION SOURCE ---
% This version GENERATES a SYNTHETIC chirp waveform (10-200 kHz) internally.
%
% Date: 2025-05-29
%

clearvars; clc; close all;

% --- Initialize Field II ---
field_init(-1);

% --- Transducer Parameters (Temporal) ---
fc_nominal = 1.0e5; % Nominal Center Frequency [Hz] (Used for lambda reference)
c = 1540;           % Speed of Sound [m/s]
lambda = c/fc_nominal; % Wavelength [m] (based on nominal fc)
fs = 2e6;           % Simulation Sampling Frequency [Hz]. High fixed value for chirp.
fprintf('Simulation sampling frequency fs = %g MHz\n', fs/1e6);

% --- Transducer Parameters (Spatial - for 3-pMUT Array) ---
pMUT_width = 1.5*lambda;    % Width of each pMUT element
pMUT_height = pMUT_width;   % Height of each pMUT element
kerf = 0.5*lambda;          % Kerf between pMUT elements (if they were in a dense grid)
d_spacing = 5*lambda;       % Characteristic spacing for the triangular pMUT array

% Define desired positions for a 3-element triangular array
pos1 = [d_spacing,0,0];
pos2 = [d_spacing*cos(2*pi/3), d_spacing*sin(2*pi/3),0];
pos3 = [d_spacing*cos(4*pi/3), d_spacing*sin(4*pi/3),0];
desired_positions = [pos1; pos2; pos3];
num_active_intended = size(desired_positions, 1);
fprintf('Intending to use %d pMUTs in a triangular configuration.\n', num_active_intended);

% Define a Field II grid large enough to encompass these desired positions
array_span_x = 2.5*d_spacing; % Span of the encompassing grid
array_span_y = 2.5*d_spacing;
element_width_grid = pMUT_width;  % Grid element size matches pMUT size
element_height_grid = pMUT_height;

num_x_grid = ceil(array_span_x / (element_width_grid + kerf));
num_y_grid = ceil(array_span_y / (element_height_grid + kerf));
if mod(num_x_grid, 2) == 0; num_x_grid = num_x_grid + 1; end % Ensure odd for centering if needed
if mod(num_y_grid, 2) == 0; num_y_grid = num_y_grid + 1; end

% --- Voxels (Imaging Grid) ---
grid_width = 2*lambda;         % Lateral width of the imaging region
grid_depth_start = 12*lambda;   % Starting depth of the imaging region
grid_depth_end = 40*lambda;     % Ending depth of the imaging region
grid_step = lambda/2;           % Pixel size (resolution)

x_coords_img = -grid_width/2 : grid_step : grid_width/2;
z_coords_img = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
Y_mesh = zeros(size(X_mesh)); % For 2D imaging in x-z plane
N_pixels = numel(X_mesh); % Total number of pixels in the image

fprintf('Wavelength lambda (nominal) = %g mm\n', lambda*1e3);
fprintf('Imaging Grid: %d pixels (%d axial x %d lateral).\n', N_pixels, length(z_coords_img), length(x_coords_img));

% --- Set Field Parameters ---
set_field('fs', fs);
set_field('c', c);

% --- Define Impulse Response and Global Excitation ---
% Using a delta function for impulse response to model a broadband pMUT
impulseResponse = [1];
fprintf('Set impulse response to a delta function (broadband transducer model).\n');
excitationPulse = 1; % Global excitation (largely overridden by ele_waveform)

% --- Define pMUT Aperture (Tx and Rx are the same) ---
% Map desired pMUT positions to the Field II grid
physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
element_no_grid_map = 0;
center_offset_x = (num_x_grid - 1)/2 * (element_width_grid + kerf);
center_offset_y = (num_y_grid - 1)/2 * (element_height_grid + kerf);
for iy = 1:num_y_grid
    y_pos_el = (iy-1)*(element_height_grid + kerf) - center_offset_y;
    for ix = 1:num_x_grid
        x_pos_el = (ix-1)*(element_width_grid + kerf) - center_offset_x;
        element_no_grid_map = element_no_grid_map + 1;
        physical_element_centers(element_no_grid_map, :) = [x_pos_el, y_pos_el, 0];
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
    warning('Could not map to %d unique elements (found %d). Check pMUT positions and grid definition.', num_active_intended, num_active);
    if num_active == 0; error('No active elements found!'); end
end
fprintf('Successfully mapped %d unique active pMUTs.\n', num_active);

enabled_matrix = zeros(num_y_grid, num_x_grid);
[row_indices, col_indices] = ind2sub([num_y_grid, num_x_grid], active_indices_linear);
for i = 1:num_active
    enabled_matrix(row_indices(i), col_indices(i)) = 1;
end

% Create the pMUT aperture in Field II
% For pMUTs, Tx and Rx are typically the same aperture.
nominal_focus_point = [0 0 100]; % Far focus, less critical with ele_waveform
sub_x_divs = 1; % No. of mathematical subdivisions of elements
sub_y_divs = 1;
pMUT_Aperture = xdc_2d_array(num_x_grid, num_y_grid, element_width_grid, element_height_grid, kerf, kerf, enabled_matrix, sub_x_divs, sub_y_divs, nominal_focus_point);

xdc_impulse(pMUT_Aperture, impulseResponse);
xdc_excitation(pMUT_Aperture, excitationPulse);

% --- Plot pMUT Array Geometry (Figure 1) ---
figure(1); clf;
[rows_active_plot, cols_active_plot] = find(enabled_matrix);
active_centers_plot = zeros(num_active, 3);
for i = 1:num_active
    iy_plot = rows_active_plot(i); ix_plot = cols_active_plot(i);
    y_pos_plot = (iy_plot-1)*(element_height_grid + kerf) - center_offset_y;
    x_pos_plot = (ix_plot-1)*(element_width_grid + kerf) - center_offset_x;
    active_centers_plot(i, :) = [x_pos_plot, y_pos_plot, 0];
end
plot3(physical_element_centers(:,1)*1e3, physical_element_centers(:,2)*1e3, physical_element_centers(:,3)*1e3, 'k.', 'MarkerSize', 1, 'DisplayName', 'Grid Elements');
hold on;
plot3(active_centers_plot(:,1)*1e3, active_centers_plot(:,2)*1e3, active_centers_plot(:,3)*1e3, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Active pMUTs');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
title(sprintf('pMUT Array Geometry (%d Active Elements)', num_active));
axis equal; grid on; view(2); legend('Location', 'best'); set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig1_pMUT_Geometry.png');

%% Generate Coded Acquisitions and Measure H components
R_acquisitions = 10; % Number of coded acquisitions (was 20, reduced for speed)
fprintf('\n--- Generating H Matrix Components from %d Coded Acquisitions ---\n', R_acquisitions);

hydrophone_positions_img = [X_mesh(:), Y_mesh(:), Z_mesh(:)]; % Imaging points

% --- Define Synthetic Chirp for Excitation ---
f_start_chirp = 10e3;  % 10 kHz
f_end_chirp = 200e3; % 200 kHz
chirp_duration = 0.2e-3; % 0.2 ms duration (was 1.5ms, reduced for speed -> smaller K_global)
fprintf('Synthetic chirp: %g-%g kHz, duration %g ms.\n', f_start_chirp/1e3, f_end_chirp/1e3, chirp_duration*1e3);

t_chirp_vec = 0 : 1/fs : chirp_duration;
if isempty(t_chirp_vec) || length(t_chirp_vec) < 2 % Ensure t_chirp_vec is valid for chirp
    t_chirp_vec = [0 1/fs]; % Minimal valid vector if duration is too short for fs
    warning('Chirp duration very short for fs, using minimal time vector for chirp.');
end
synth_chirp_base = chirp(t_chirp_vec, f_start_chirp, t_chirp_vec(end), f_end_chirp, 'linear');
synth_chirp_windowed = synth_chirp_base .* tukeywin(length(t_chirp_vec), 0.25)'; % Apply window
base_chirp_to_modulate = synth_chirp_windowed(:)'; % Ensure row vector

% Store data from each acquisition
all_hhp_data = cell(R_acquisitions, 1);
all_start_times = zeros(R_acquisitions, 1);
all_K_values = zeros(R_acquisitions, 1); % Number of time samples for each hhp_r

rng('default'); % For reproducible random modulations

% --- Plot Example Transmitted Field (Figure 2 adaptation) ---
% (This is a simplified version, showing field from one pMUT with the base chirp)
figure(2); clf;
if num_active > 0
    ele_waveform(pMUT_Aperture, 1, base_chirp_to_modulate); % Excite only the first pMUT
    for p_idx_other = 2:num_active % Silence other pMUTs for this specific plot
        ele_waveform(pMUT_Aperture, p_idx_other, 0);
    end
    % Define a line of points in front of the first pMUT for field calculation
    first_pmut_pos = active_centers_plot(1,:);
    field_calc_x = first_pmut_pos(1); % Keep x constant
    field_calc_y = first_pmut_pos(2); % Keep y constant
    field_calc_z = (grid_depth_start : lambda/4 : grid_depth_end);
    field_calc_points = [repmat(field_calc_x, length(field_calc_z),1), ...
                         repmat(field_calc_y, length(field_calc_z),1), ...
                         field_calc_z(:)];
    if ~isempty(field_calc_points)
        [hp_example, start_hp_example] = calc_hp(pMUT_Aperture, field_calc_points);
        hp_env_example = abs(hilbert(hp_example));
        imagesc(t_chirp_vec*1e6, field_calc_z*1e3, hp_env_example'); % Transpose for correct orientation
        xlabel('Time (us)'); ylabel('Depth (mm)');
        title('Example Transmitted Field (Envelope from 1st pMUT)');
        colorbar;
    else
        title('Example Transmitted Field (Skipped - No field points)');
    end
else
    title('Example Transmitted Field (Skipped - No active pMUTs)');
end
sgtitle('Figure 2: Example Transmitted Field');
set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig2_ExampleTxField.png');
% Reset waveforms for actual H generation
for p_idx = 1:num_active
    ele_waveform(pMUT_Aperture, p_idx, 0); % Clear previous waveforms
end


% Loop for each coded acquisition
for r_acq = 1:R_acquisitions
    fprintf('Calculating Acquisition %d/%d...\n', r_acq, R_acquisitions);
    % Apply unique, phase-shifted chirps to each pMUT
    for p_idx = 1:num_active % Iterate through the active pMUTs
        current_rand_phase = 2 * pi * rand(); % Random phase for this pMUT, this acquisition
        
        analytic_signal = hilbert(base_chirp_to_modulate);
        phase_shifted_analytic = analytic_signal * exp(1i * current_rand_phase);
        current_pmut_waveform = real(phase_shifted_analytic);
        
        ele_waveform(pMUT_Aperture, p_idx, current_pmut_waveform);
    end
    
    % Calculate pulse-echo field for this coded transmission
    [hhp_r, start_time_r] = calc_hhp(pMUT_Aperture, pMUT_Aperture, hydrophone_positions_img);
    
    all_hhp_data{r_acq} = hhp_r;
    all_start_times(r_acq) = start_time_r;
    all_K_values(r_acq) = size(hhp_r, 1);
    
    if ~isempty(hhp_r)
        fprintf('    Max abs value in hhp_r: %g (K_r = %d)\n', max(abs(hhp_r(:))), size(hhp_r,1));
    else
        fprintf('    WARNING: hhp_r is empty for acquisition %d!\n', r_acq);
    end
end

%% Assemble Full H Matrix with Interpolation and Alignment
fprintf('\n--- Assembling Full H Matrix ---\n');
all_end_times = zeros(R_acquisitions, 1);
for r_idx = 1:R_acquisitions
    if all_K_values(r_idx) > 0
        all_end_times(r_idx) = all_start_times(r_idx) + (all_K_values(r_idx) - 1)/fs;
    else
        all_end_times(r_idx) = all_start_times(r_idx); % Handle case of empty hhp_r
    end
end

min_global_start_time = min(all_start_times);
max_global_end_time = max(all_end_times);

if min_global_start_time >= max_global_end_time % Check if time window is valid
     if any(all_K_values > 0) % If there was some data
        max_global_end_time = min_global_start_time + (max(all_K_values(all_K_values>0))-1)/fs; % Estimate based on longest K
        if min_global_start_time >= max_global_end_time % Still invalid
            max_global_end_time = min_global_start_time + 100/fs; % Fallback: 100 samples duration
            warning('Global time window for H assembly was invalid, using fallback.');
        end
     else % No data at all
        min_global_start_time = 0; max_global_end_time = 100/fs; % Arbitrary small window
        warning('No data in any acquisition, H will be all zeros.');
     end
end


t_common_axis = min_global_start_time : 1/fs : max_global_end_time;
K_global_per_acq = length(t_common_axis);

fprintf('Global Time Window for H assembly: Start = %g s, End = %g s, K_global_per_acq = %d samples.\n', ...
        min_global_start_time, max_global_end_time, K_global_per_acq);

H_assembled = zeros(K_global_per_acq * R_acquisitions, N_pixels);
current_row_offset = 0;
for r_acq = 1:R_acquisitions
    hhp_current = all_hhp_data{r_acq};
    start_time_current = all_start_times(r_acq);
    K_current = all_K_values(r_acq);
    
    if K_current == 0 || isempty(hhp_current)
        fprintf('  Skipping interpolation for acquisition %d (K=0 or empty hhp).\n', r_acq);
        current_row_offset = current_row_offset + K_global_per_acq;
        continue;
    end
    
    t_current_acq_axis = start_time_current + (0:(K_current-1))/fs;
    
    hhp_aligned_r = zeros(K_global_per_acq, N_pixels);
    for px_col = 1:N_pixels
        if ~isempty(hhp_current) && size(hhp_current,2) >= px_col
             % Ensure t_current_acq_axis is monotonically increasing for interp1
            if length(t_current_acq_axis) > 1 && issorted(t_current_acq_axis)
                hhp_aligned_r(:, px_col) = interp1(t_current_acq_axis, hhp_current(:,px_col), t_common_axis, 'linear', 0);
            elseif length(t_current_acq_axis) == 1 % Single point, replicate if it falls in t_common_axis
                 idx_match = find(abs(t_common_axis - t_current_acq_axis) < (0.5/fs) , 1);
                 if ~isempty(idx_match)
                    hhp_aligned_r(idx_match, px_col) = hhp_current(1,px_col);
                 end
            else
                 % t_current_acq_axis is not sorted or invalid, fill with zeros
            end
        end
    end
    
    H_assembled(current_row_offset + (1:K_global_per_acq), :) = hhp_aligned_r;
    current_row_offset = current_row_offset + K_global_per_acq;
end

H = H_assembled;
M = size(H, 1); % Total number of measurements (rows in H)
N = N_pixels;   % Total number of pixels (columns in H)

fprintf('Final Assembled H Matrix: %d rows (M) x %d columns (N).\n', M, N);
if M == 0; error('H matrix has zero rows. Cannot proceed.'); end

% --- Plot a few columns of H (Figure 3 adaptation) ---
figure(3); clf;
hold on;
num_cols_to_plot_H = min(N, 5);
indices_to_plot_H = round(linspace(1, N, num_cols_to_plot_H));
% Time axis for the entire assembled H matrix (concatenated acquisitions)
t_axis_H_plot_full_concat = (0:(M-1))/fs * 1e6; % Time in microseconds

for n_idx_plot = 1:length(indices_to_plot_H)
    col_idx = indices_to_plot_H(n_idx_plot);
    plot(t_axis_H_plot_full_concat, H(:, col_idx), 'DisplayName', sprintf('H column for Pixel %d', col_idx));
end
hold off;
xlabel('Overall Row Index (Effective Concatenated Time - us)');
ylabel('Amplitude');
title('Figure 3: Columns of Assembled H Matrix');
axis tight; grid on; legend('Location', 'best'); set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig3_H_Columns.png');

%% Generate a Scene to be Imaged (v)
fprintf('\n--- Generating Scene (v) ---\n');
scene_xmin = min(x_coords_img); scene_xmax = max(x_coords_img);
scene_zmin = min(z_coords_img); scene_zmax = max(z_coords_img);
scene_z_range = scene_zmax - scene_zmin;

% Define targets (positions and amplitudes)
targets.x = [scene_xmin*0.5, 0, scene_xmax*0.5, scene_xmin*0.2];
targets.y = zeros(size(targets.x)); % 2D imaging
targets.z = [scene_zmin + 0.2*scene_z_range, ...
             scene_zmin + 0.5*scene_z_range, ...
             scene_zmin + 0.8*scene_z_range, ...
             scene_zmin + 0.4*scene_z_range];
amplitudes = [1, 0.9, 1.1, 0.7]'; % Ensure column vector

scene_matrix = zeros(length(z_coords_img), length(x_coords_img));
for i = 1:length(targets.x)
    [~, ix_scene] = min(abs(x_coords_img - targets.x(i)));
    [~, iz_scene] = min(abs(z_coords_img - targets.z(i)));
    if ix_scene > 0 && ix_scene <= length(x_coords_img) && iz_scene > 0 && iz_scene <= length(z_coords_img)
        scene_matrix(iz_scene, ix_scene) = amplitudes(i);
    else
        warning('Target %d at (%.2f mm, %.2f mm) is outside the defined grid.', i, targets.x(i)*1e3, targets.z(i)*1e3);
    end
end
v_true_vector = scene_matrix(:); % Ground truth image vector (N x 1)
fprintf('Scene generated with %d targets.\n', length(targets.x));

% --- Plot True Scene (Figure 4) ---
figure(4); clf;
subplot(1, 2, 1);
imagesc(x_coords_img*1e3, z_coords_img*1e3, scene_matrix);
axis image; colormap gray; colorbar;
xlabel('x (mm)'); ylabel('z (mm)'); title('True Image v (matrix form)');
set(gca, 'YDir','normal');

subplot(1, 2, 2);
stem(v_true_vector);
axis square tight;
title('True Image v (vector form)'); xlabel('Pixel Index'); ylabel('Amplitude');
sgtitle('Figure 4: True Image (Ground Truth)');
set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig4_TrueImage.png');

%% Image Formation Model: u = H*v + n
fprintf('\n--- Simulating Measurement Acquisition (u = Hv + n) ---\n');

if isempty(H) || isempty(v_true_vector)
    error('H or v_true_vector is empty. Cannot proceed.');
end
if size(H,2) ~= length(v_true_vector)
    error('Dimension mismatch: H columns (%d) must match length of v (%d).', size(H,2), length(v_true_vector));
end

Hv_signal = H * v_true_vector; % Ideal signal (M x 1)

% Add Gaussian noise
rng('shuffle'); % Use 'shuffle' for different noise each run
electronic_SNR_db = 30; % Define desired SNR in dB
electronic_SNR = 10^(electronic_SNR_db/10);
signal_power_est = mean(Hv_signal(:).^2);

if isnan(signal_power_est) || signal_power_est < 1e-20 % Threshold for very weak signal
    warning('Signal power is effectively zero or NaN (%g). Using small default noise sigma.', signal_power_est);
    noise_sigma = 1e-7; % Arbitrary small noise if signal is too weak
    actual_SNR_db = -Inf;
else
    noise_variance = signal_power_est / electronic_SNR;
    noise_sigma = sqrt(noise_variance);
    actual_SNR_db = 10*log10(signal_power_est / noise_variance); % Actual SNR
end

n_noise_vec = noise_sigma * randn(size(Hv_signal)); % Noise vector (M x 1)
u_measured_signal = Hv_signal + n_noise_vec; % Measured signal (M x 1)

fprintf('Added Gaussian noise. Target SNR: %.1f dB. Actual SNR: %.1f dB. Noise sigma: %g\n', ...
        electronic_SNR_db, actual_SNR_db, noise_sigma);

% --- Plot Measurements (Figure 6 - similar to original Ensemble.m) ---
figure(6); clf;
t_axis_measurements_plot = (0:(M-1))/fs * 1e6; % Time in microseconds
plot_limit_samples = min(M, round(fs*max_global_end_time*1.1)); % Plot slightly beyond max signal end
if plot_limit_samples == 0 && M > 0; plot_limit_samples = M; end % ensure some plotting if M>0

subplot(1, 3, 1);
if ~isempty(Hv_signal) && plot_limit_samples > 0; plot(t_axis_measurements_plot(1:plot_limit_samples), Hv_signal(1:plot_limit_samples)); else; plot([]); end
axis tight; grid on; xlabel('Time (µs)'); title('Measurements (Hv - Ideal)'); ylabel('Amplitude');
yl_hv_plot = ylim;

subplot(1, 3, 2);
if ~isempty(n_noise_vec) && plot_limit_samples > 0; plot(t_axis_measurements_plot(1:plot_limit_samples), n_noise_vec(1:plot_limit_samples)); else; plot([]); end
axis tight; grid on; xlabel('Time (µs)'); title(sprintf('Noise (SNR \\approx %.1f dB)', actual_SNR_db)); ylabel('Amplitude');
if ~isempty(yl_hv_plot) && all(isfinite(yl_hv_plot)); ylim(yl_hv_plot); end

subplot(1, 3, 3);
if ~isempty(u_measured_signal) && plot_limit_samples > 0; plot(t_axis_measurements_plot(1:plot_limit_samples), u_measured_signal(1:plot_limit_samples)); else; plot([]); end
axis tight; grid on; xlabel('Time (µs)'); title('Measurements + Noise (u)'); ylabel('Amplitude');
if ~isempty(yl_hv_plot) && all(isfinite(yl_hv_plot)); ylim(yl_hv_plot); end

sgtitle(sprintf('Figure 6: Measurement Signals (Compression = %.2f)', N/M));
set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig6_MeasuredSignals.png');

%% Reconstruction Algorithms
fprintf('\n--- Starting Reconstruction ---\n');

% Define common variables for reconstruction, similar to Ensemble.m
A_matrix = H;                     % System matrix (M x N)
At_matrix = transpose(A_matrix);  % Transposed system matrix (N x M)
b_vector = u_measured_signal;     % Measurement vector (M x 1)
I_true_matrix = scene_matrix;     % True image (for comparison) (Nz x Nx)
v_true_vec_norm = v_true_vector;  % For PSNR calculation
if max(abs(v_true_vec_norm)) > 0
    v_true_vec_norm = v_true_vec_norm ./ max(abs(v_true_vec_norm));
end

imageResolution = size(I_true_matrix); % For ADMM reshaping [Nz, Nx]

% Function handles for matrix operations
Afun = @(x_vec) A_matrix*x_vec;         % x_vec is N x 1 image vector
Atfun = @(y_vec) At_matrix*y_vec;       % y_vec is M x 1 measurement vector

% --- Least Norm solution with Predetermined Conjugate Gradient (PCG) ---
fprintf('  Calculating Least Norm (PCG) Solution...\n');
tic;
% AAtfun for PCG: A*(At*y_intermediate_cg_vec)
% Input y_intermediate_cg_vec is M x 1, output is M x 1
AAtfun_pcg  = @(y_intermediate_cg_vec) Afun(Atfun(y_intermediate_cg_vec));

maxItersCG_main = 200; % Max iterations for PCG (Original Ensemble.m used 1000)
tolCG_main = noise_sigma; % Tolerance based on noise estimation

if M==0 || N==0 || isempty(b_vector) || norm(b_vector)==0
    fprintf('    Skipping PCG (Zero Dim M or N, or empty/zero b_vector).\n');
    x_pcg_img_vec = zeros(N,1); iter_cg_main = 0; flag_cg_main = -1;
else
    % PCG solves (A A^T)y = b, then x_estimate = A^T y
    [y_cg_solution_M_by_1, flag_cg_main, ~, iter_cg_main] = pcg(AAtfun_pcg, b_vector(:), tolCG_main, maxItersCG_main);
    if flag_cg_main ~= 0; warning('PCG did not converge (Flag=%d, Iters=%d). Result might be inaccurate.', flag_cg_main, iter_cg_main); end
    x_pcg_img_vec = Atfun(y_cg_solution_M_by_1); % Image estimate (N x 1)
end

% Post-process and normalize PCG result
x_pcg_img_vec_processed = real(x_pcg_img_vec);
min_val_pcg = min(x_pcg_img_vec_processed(:));
max_val_pcg = max(x_pcg_img_vec_processed(:));
if (max_val_pcg - min_val_pcg) > eps
    x_pcg_norm_vec = (x_pcg_img_vec_processed - min_val_pcg) / (max_val_pcg - min_val_pcg);
else
    x_pcg_norm_vec = zeros(N,1); % Avoid NaN if all values are same
end
x_pcg_norm_vec(isnan(x_pcg_norm_vec)) = 0;
x_pcg2D = reshape(x_pcg_norm_vec, imageResolution);

MSE_pcg = mean( (x_pcg_norm_vec - v_true_vec_norm).^2 );
PSNR_pcg = 10*log10(1/MSE_pcg);
runtime_pcg = toc;
fprintf('  PCG: MSE = %g, PSNR = %.2f dB, Iters = %d, Time = %.2f s\n', MSE_pcg, PSNR_pcg, iter_cg_main, runtime_pcg);

% --- Least Norm solution with Moore-Penrose Pseudo-inverse ---
fprintf('  Calculating Least Norm (Pseudo-inverse) Solution...\n');
tic;
tolPinv = 1e-6; % Tolerance for pinv
x_pinv_img_vec = zeros(N, 1); % Initialize
pinv_method_desc = 'Skipped';

if M > 0 && N > 0 && ~isempty(A_matrix) && ~isempty(b_vector)
    try
        if M < N % Underdetermined (more pixels than measurements per acquisition time point)
            AAt_pinv_calc = A_matrix*At_matrix; % M x M matrix
            if rank(AAt_pinv_calc) < min(size(AAt_pinv_calc))
                fprintf('    Warning: A*At is rank deficient for Pinv. Rank = %d.\n', rank(AAt_pinv_calc));
            end
            AAt_inv = pinv(AAt_pinv_calc, tolPinv); % M x M
            x_pinv_img_vec = At_matrix * (AAt_inv * b_vector(:)); % N x 1
            pinv_method_desc = 'At*pinv(A*At)*b';
        else % Overdetermined or M=N
            A_pinv_calc = pinv(A_matrix, tolPinv); % N x M
            x_pinv_img_vec = A_pinv_calc * b_vector(:); % N x 1
            pinv_method_desc = 'pinv(A)*b';
        end
        fprintf('    Used Pinv method: %s\n', pinv_method_desc);
    catch ME
        fprintf('    Pinv calculation failed: %s\n', ME.message);
        pinv_method_desc = 'Failed';
    end
else
     fprintf('    Skipping Pinv (Zero Dim M or N, or empty A/b_vector).\n');
end

% Post-process and normalize Pinv result
x_pinv_img_vec_processed = real(x_pinv_img_vec);
min_val_pinv = min(x_pinv_img_vec_processed(:));
max_val_pinv = max(x_pinv_img_vec_processed(:));
if (max_val_pinv - min_val_pinv) > eps
    x_pinv_norm_vec = (x_pinv_img_vec_processed - min_val_pinv) / (max_val_pinv - min_val_pinv);
else
    x_pinv_norm_vec = zeros(N,1);
end
x_pinv_norm_vec(isnan(x_pinv_norm_vec))=0;
x_pinv2D = reshape(x_pinv_norm_vec, imageResolution);

MSE_pinv = mean( (x_pinv_norm_vec - v_true_vec_norm).^2 );
PSNR_pinv = 10*log10(1/MSE_pinv);
runtime_pinv = toc;
fprintf('  Pinv (%s): MSE = %g, PSNR = %.2f dB, Time = %.2f s\n', pinv_method_desc, MSE_pinv, PSNR_pinv, runtime_pinv);

% --- Plot Least Norm Results (Figure 7) ---
figure(7); clf;
subplot(1, 3, 1);
imagesc(x_coords_img*1e3, z_coords_img*1e3, I_true_matrix);
axis image; colormap gray; colorbar; set(gca, 'YDir','normal');
title('Ground Truth'); xlabel('x (mm)'); ylabel('z (mm)');

subplot(1, 3, 2);
imagesc(x_coords_img*1e3, z_coords_img*1e3, x_pcg2D);
axis image; colormap gray; colorbar; set(gca, 'YDir','normal');
title(sprintf('Least Norm (PCG)\nIters=%d, PSNR=%.2fdB', iter_cg_main, PSNR_pcg));
xlabel('x (mm)'); ylabel('z (mm)');

subplot(1, 3, 3);
imagesc(x_coords_img*1e3, z_coords_img*1e3, x_pinv2D);
axis image; colormap gray; colorbar; set(gca, 'YDir','normal');
title(sprintf('Least Norm (Pinv %s)\nPSNR=%.2fdB', pinv_method_desc, PSNR_pinv));
xlabel('x (mm)'); ylabel('z (mm)');

sgtitle(sprintf('Figure 7: Least Norm Reconstructions (R_{acq}=%d, SNR~%.1fdB)', R_acquisitions, actual_SNR_db));
set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig7_LeastNormRecon.png');

%% ADMM with TV regularization (similar to Ensemble.m)
fprintf('  Calculating ADMM-TV Solution...\n');

% Normalize system matrix A_matrix for ADMM stability
H_norm_factor = max(abs(A_matrix(:)));
if H_norm_factor < eps; H_norm_factor = 1; end % Avoid division by zero
A_admm = A_matrix ./ H_norm_factor;
At_admm = transpose(A_admm);

% Scaled measurement vector and noise sigma for normalized system
b_admm_vec = b_vector(:) / H_norm_factor;
noise_sigma_admm = noise_sigma / H_norm_factor; % Scale noise_sigma too

% ADMM function handles (operating on images or vectors as appropriate)
Afun_admm = @(x_img_vec) A_admm * x_img_vec; % Input: N x 1 image vector; Output: M x 1
Atfun_admm_img = @(y_vec) reshape(At_admm * y_vec, imageResolution); % Input: M x 1; Output: Nz x Nx image
AtAfun_admm_img = @(x_img) Atfun_admm_img(Afun_admm(x_img(:))); % Input: Nz x Nx image; Output: Nz x Nx image

% Define difference operators for TV (from helper function at the end)
[Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
% opDx: Input image (Nz x Nx), Output [N_pixels x 2] (Dx*x_vec, Dy*x_vec)
opDx_tv = @(x_img_tv) [Dx_sparse*x_img_tv(:), Dy_sparse*x_img_tv(:)];
% opDtx: Input [N_pixels x 2], Output image (Nz x Nx)
opDtx_tv = @(v_grad_vecs_tv) reshape(Dx_sparse'*v_grad_vecs_tv(:,1) + Dy_sparse'*v_grad_vecs_tv(:,2), imageResolution);
% opDtDx: Input image (Nz x Nx), Output image (Nz x Nx)
opDtDx_tv = @(x_img_tv) opDtx_tv(opDx_tv(x_img_tv));

% ADMM parameters (from original Ensemble.m where applicable)
numItersADMM    = 25;  % Number of ADMM iterations
rho_admm        = 10;  % ADMM penalty parameter
lambda_tv_reg   = 0.01;% TV regularization weight

% Initialize ADMM variables
x_admm_img_iter = zeros(imageResolution); % Current image estimate (Nz x Nx)
z_admm_grad_iter = zeros([prod(imageResolution) 2]); % Auxiliary for gradients (N_pixels x 2)
u_admm_dual_iter = zeros([prod(imageResolution) 2]); % Dual variable for gradients (N_pixels x 2)

% Function handle for PCG within ADMM x-update (original style)
% Hfun_admm_pcg: (AtA + rho*DtD) * x_vec
Hfun_pcg_admm = @(x_img_vec_pcg) ...
    reshape( AtAfun_admm_img(reshape(x_img_vec_pcg,imageResolution)) + ...
             rho_admm.*opDtDx_tv(reshape(x_img_vec_pcg,imageResolution)), ...
             [prod(imageResolution) 1] );

% Store results per iteration
PSNR_admm_iters = zeros([numItersADMM 1]);
residuals_admm_iters = zeros([numItersADMM 1]); % For objective function value

% --- ADMM Loop (Figure 8 plotting) ---
figure(8); clf;
tic;
for k_admm = 1:numItersADMM
    % x update (image estimate)
    % Solve (AtA + rho*DtD)x = At*b_admm + rho*Dt(z-u)
    v_for_x_update = z_admm_grad_iter - u_admm_dual_iter;
    bb_for_x_update = Atfun_admm_img(b_admm_vec) + rho_admm*opDtx_tv(v_for_x_update);
    
    maxItersCG_admm = 25; % Max iterations for inner PCG
    % Use scaled noise_sigma_admm as tolerance for inner PCG (like original Ensemble.m)
    tolCG_admm = noise_sigma_admm; 
    
    if norm(bb_for_x_update(:)) == 0 % If RHS is zero, solution is zero
        x_admm_img_vec_new = zeros(prod(imageResolution),1);
        flag_pcg_admm_inner = 0; iter_pcg_admm_inner = 0;
    else
        [x_admm_img_vec_new, flag_pcg_admm_inner, ~, iter_pcg_admm_inner] = ...
            pcg(Hfun_pcg_admm, bb_for_x_update(:), tolCG_admm, maxItersCG_admm, [], [], x_admm_img_iter(:));
    end
    
    if flag_pcg_admm_inner ~= 0 && iter_pcg_admm_inner > 0 % Don't warn if iter is 0 (e.g. tol met instantly)
        fprintf('    ADMM Iter %d: Inner PCG Flag=%d, Iters=%d.\n', ...
                k_admm, flag_pcg_admm_inner, iter_pcg_admm_inner);
    end
    x_admm_img_iter = reshape(x_admm_img_vec_new, imageResolution);
    
    % z update (gradient shrinkage - vectorial for isotropic TV)
    kappa_tv_shrink = lambda_tv_reg / rho_admm;
    v_for_z_update = opDx_tv(x_admm_img_iter) + u_admm_dual_iter; % N_pixels x 2
    
    v_grad_norm_z = sqrt(sum(v_for_z_update.^2, 2)); % N_pixels x 1
    v_grad_norm_z(v_grad_norm_z < eps) = eps; % Avoid division by zero
    
    shrinkage_factor_z = max(0, 1 - kappa_tv_shrink ./ v_grad_norm_z); % N_pixels x 1
    z_admm_grad_iter = v_for_z_update .* shrinkage_factor_z; % N_pixels x 2
    
    % u update (dual variable)
    u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
    
    % Metrics for this iteration
    % Scale x_admm_img_iter for display/PSNR to [0, 1] range
    x_admm_scaled_vec = real(x_admm_img_iter(:));
    min_x_admm_iter = min(x_admm_scaled_vec); max_x_admm_iter = max(x_admm_scaled_vec);
    if (max_x_admm_iter - min_x_admm_iter) > eps
        x_admm_scaled_vec = (x_admm_scaled_vec - min_x_admm_iter) / (max_x_admm_iter - min_x_admm_iter);
    else
        x_admm_scaled_vec = zeros(size(x_admm_scaled_vec));
    end
    x_admm_scaled_vec(isnan(x_admm_scaled_vec)) = 0;

    MSE_current_admm = mean( (x_admm_scaled_vec - v_true_vec_norm).^2 );
    PSNR_admm_iters(k_admm) = 10*log10(1/MSE_current_admm);
    
    % Objective function: 0.5*||Ax-b||^2 + lambda*||Dx||_1 (using isotropic TV norm)
    r1_fidelity = b_admm_vec - Afun_admm(x_admm_img_iter(:));
    r2_tv_opDx_val = opDx_tv(x_admm_img_iter); % N_pixels x 2
    tv_norm_val = sum(sqrt(sum(r2_tv_opDx_val.^2, 2))); % Sum of L2 norms of gradient vectors
    residuals_admm_iters(k_admm) = 0.5*sum(r1_fidelity(:).^2) + lambda_tv_reg * tv_norm_val;
    
    runtime_ADMM_iter_curr = toc;
    
    % Plotting update (similar to Ensemble.m)
    subplot(2,3,[1 4]);
    imagesc(x_coords_img*1e3, z_coords_img*1e3, I_true_matrix); axis image; colormap(gca, gray); colorbar;
    set(gca,'YDir','normal'); title('Target Image'); xlabel('x (mm)'); ylabel('z (mm)');
    
    subplot(2,3,[2 5]);
    imagesc(x_coords_img*1e3, z_coords_img*1e3, reshape(x_admm_scaled_vec, imageResolution)); axis image; colormap(gca, gray); colorbar;
    set(gca,'YDir','normal'); title(sprintf('\\lambda=%.2e, \\rho=%.1f\nPSNR=%.2fdB, Iter %d', lambda_tv_reg, rho_admm, PSNR_admm_iters(k_admm), k_admm));
    xlabel('x (mm)'); ylabel('z (mm)');
    
    subplot(2,3,3);
    plot(1:k_admm, PSNR_admm_iters(1:k_admm), 'LineWidth', 2, 'color', [1 0 1]);
    title('PSNR vs Iteration'); xlabel('Iteration'); ylabel('PSNR (dB)'); grid on; axis tight;
    if k_admm > 1; current_psnr_min = min(PSNR_admm_iters(1:k_admm)); current_psnr_max = max(PSNR_admm_iters(1:k_admm)); if current_psnr_min == current_psnr_max; current_psnr_max = current_psnr_min+1; end; if isfinite(current_psnr_min) && isfinite(current_psnr_max); ylim([current_psnr_min-1, current_psnr_max+1]); end; end

    subplot(2,3,6);
    plot(1:k_admm, log10(residuals_admm_iters(1:k_admm)), 'LineWidth', 2);
    title('log10(Objective) vs Iteration'); xlabel('Iteration'); ylabel('log10(Value)'); grid on; axis tight;
    if k_admm > 1; current_res_min = min(log10(residuals_admm_iters(1:k_admm))); current_res_max = max(log10(residuals_admm_iters(1:k_admm))); if current_res_min == current_res_max; current_res_max = current_res_min+1; end; if isfinite(current_res_min) && isfinite(current_res_max); ylim([current_res_min-0.5, current_res_max+0.5]); end; end
    
    drawnow;
    if mod(k_admm,5)==0 || k_admm == numItersADMM
      fprintf('  ADMM Iter %d/%d: PSNR = %.2f dB, log10(Obj) = %.3f, Runtime = %.2f s\n', ...
              k_admm, numItersADMM, PSNR_admm_iters(k_admm), log10(residuals_admm_iters(k_admm)), runtime_ADMM_iter_curr);
    end
end % End ADMM iterations loop
runtime_ADMM_total = toc;

fprintf('  ADMM-TV complete: Final PSNR = %.2f dB, Total Time = %.2fs\n', PSNR_admm_iters(end), runtime_ADMM_total);
sgtitle(sprintf('Figure 8: ADMM TV (R_{acq}=%d, SNR~%.1fdB)', R_acquisitions, actual_SNR_db));
set(gcf, 'Color', 'w');
saveas(gcf, 'SynthChirp_Fig8_ADMM_TV.png');

%% End Field II Simulation
if exist('field_end','file') == 2
    field_end;
    disp('Field II ended.');
else
    disp('field_end function not found. Field II may not have been properly initialized or path is not set.');
end

%% Helper function for TV difference operators (from Ensemble.m style)
function [Dx, Dy] = createDifferenceOperators(imageSize)
    rows = imageSize(1); % Nz
    cols = imageSize(2); % Nx
    N_img_pixels = rows * cols;

    % Difference operator in x-direction (across columns of the image matrix)
    % Corresponds to shifting by 'rows' in the vectorized image
    Dx = spdiags([-ones(N_img_pixels,1), ones(N_img_pixels,1)], [0, rows], N_img_pixels, N_img_pixels);
    % Boundary condition: Zero difference for elements in the last 'column' of the image
    % (when image is vectorized column-wise, these are the last 'rows' elements if not careful)
    % More robust: identify elements in the last column of the image matrix
    last_col_indices_mask = false(N_img_pixels,1);
    last_col_indices_mask( (cols-1)*rows + 1 : cols*rows ) = true;
    Dx(last_col_indices_mask, :) = 0; % Set rows of Dx corresponding to last image column to zero

    % Difference operator in y-direction (across rows of the image matrix)
    % Corresponds to shifting by '1' in the vectorized image
    Dy = spdiags([-ones(N_img_pixels,1), ones(N_img_pixels,1)], [0, 1], N_img_pixels, N_img_pixels);
    % Boundary condition: Zero difference for elements in the last 'row' of the image
    last_row_indices_mask = false(N_img_pixels,1);
    last_row_indices_mask(rows:rows:N_img_pixels) = true; % Elements at the end of each column
    Dy(last_row_indices_mask, :) = 0;
end
