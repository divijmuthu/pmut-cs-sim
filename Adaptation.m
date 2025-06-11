% Merged_pMUT_OriginalRecon_SimplePulse.m
%
% Combines features from Ensemble.m and Ensemble_pMUT_ChirpPhase_Workaround.m
%
% Goal: Use the pMUT array setup (triangular via 'enabled' matrix),
%       phase shifting via ele_delay.
%       Use a SIMPLE excitation pulse (like original Ensemble.m) - NO CHIRP.
%       Use the reconstruction algorithms (PCG, Pinv, ADMM) from the
%       'original' script (Ensemble.m).
%
% Version: 2025-05-04 (Simple Pulse)

clearvars; clc; close all;

% Initialize Field II
field_init(-1);

% --- Simulation Parameters ---
fc = 5.8e4;       % Carrier Frequency [Hz] - Center frequency for delays
c = 1540;       % Speed of Sound [m/s]
lambda = c/fc;  % Wavelength [m]
fBW = 0.66;      % Fractional Bandwidth (used for impulse response)
fs = fc*4;      % Sampling Frequency [Hz]

% Set Field Parameters
set_field('fs', fs); % Set Field II sampling rate
set_field('c', c);   % Set speed of sound

% --- Define Simple Excitation Pulse (Like Original Ensemble.m) ---
excitationPulse = 1; % Simple delta function excitation

% --- Define Impulse Response (Receive Filter / Effective Tx Pulse) ---
% Represents transducer bandwidth. Convolved with excitationPulse=1 gives itself.
tc = gauspuls('cutoff', fc, fBW, -6, -40);
t_ir = -tc : 1/fs : tc;
impulseResponse = gauspuls(t_ir, fc, fBW);
fprintf('Using simple delta excitation convolved with Gauspuls impulse response.\n');

% --- pMUT Element Parameters & Array Definition ---
pMUT_width = 1.5*lambda;
pMUT_height = pMUT_width;
kerf = 0.5*lambda;
kerf_x = kerf;
kerf_y = kerf;

% Define desired positions for a triangular array (3 elements)
d = 5*lambda; % Distance from center
pos1 = [d,0,0];
pos2 = [d*cos(2*pi/3), d*sin(2*pi/3),0];
pos3 = [d*cos(4*pi/3), d*sin(4*pi/3),0];
desired_positions = [pos1; pos2; pos3];
num_active_intended = size(desired_positions, 1); % Should be 3

% Define the grid bounds large enough to contain desired positions
array_span_x = 2.5*d;
array_span_y = 2.5*d;
element_width = pMUT_width;
element_height = pMUT_height;

% Calculate grid dimensions
num_x = ceil(array_span_x / (element_width + kerf_x));
num_y = ceil(array_span_y / (element_height + kerf_y));
if mod(num_x, 2) == 0; num_x = num_x + 1; end % Ensure odd number if needed
if mod(num_y, 2) == 0; num_y = num_y + 1; end
num_elements_total = num_x * num_y;
fprintf('Defining %dx%d array grid to encompass desired elements...\n', num_x, num_y);

% --- Imaging Grid ---
grid_width = 20*lambda;
grid_depth_start = 10*lambda;
grid_depth_end = 40*lambda;
grid_step = lambda/2;
x_grid = -grid_width/2 : grid_step : grid_width/2;
z_grid = grid_depth_start : grid_step : grid_depth_end;
[X, Z] = meshgrid(x_grid, z_grid);
Y = zeros(size(X));
N_pixels = numel(X);
fprintf('Imaging Grid: %d pixels...\n', N_pixels);

% --- Find Element Indices & Create 'enabled' Matrix ---
physical_element_centers = zeros(num_elements_total, 3);
element_no = 0;
center_offset_x = (num_x - 1)/2 * (element_width + kerf_x);
center_offset_y = (num_y - 1)/2 * (element_height + kerf_y);
for iy = 1:num_y
    y_pos = (iy-1)*(element_height + kerf_y) - center_offset_y;
    for ix = 1:num_x
        x_pos = (ix-1)*(element_width + kerf_x) - center_offset_x;
        element_no = element_no + 1;
        physical_element_centers(element_no, :) = [x_pos, y_pos, 0];
    end
end
active_indices_linear = zeros(num_active_intended, 1);
for i = 1:num_active_intended
    distances = sqrt(sum((physical_element_centers - desired_positions(i,:)).^2, 2));
    [~, min_idx] = min(distances);
    active_indices_linear(i) = min_idx;
end
active_indices_linear = unique(active_indices_linear); % Ensure uniqueness
num_active = length(active_indices_linear);
if num_active ~= num_active_intended
    warning('Could not map to %d unique elements (found %d). Adjust grid or desired positions if needed.', num_active_intended, num_active);
    if num_active == 0; error('No active elements found! Check grid definition and desired positions.'); end
end
fprintf('Using %d unique active element linear indices: %s\n', num_active, mat2str(active_indices_linear'));
enabled = zeros(num_y, num_x);
[row_indices, col_indices] = ind2sub([num_y, num_x], active_indices_linear);
for i = 1:num_active; enabled(row_indices(i), col_indices(i)) = 1; end
fprintf('Created enabled matrix.\n');

% --- Define the 2D Array Aperture using the 'enabled' matrix ---
focus_point = [0 0 100]; % Far focus (less critical when using delays/phasing)
sub_x = 1; % Element subdivisions for Field II
sub_y = 1;
Th = xdc_2d_array(num_x, num_y, element_width, element_height, kerf_x, kerf_y, enabled, sub_x, sub_y, focus_point);
fprintf('Created 2D array aperture handle Th using %d active elements.\n', num_active);

% --- Define Desired Phase Shifts and Calculate Time Delays ---
if num_active == 3 % Apply specific phases if we found exactly 3 elements
    target_phases_rad = [0, pi/2, pi]; % Example phases
    fprintf('Applying phases [0, pi/2, pi] to the %d active elements.\n', num_active);
elseif num_active > 0 % Default to zero phase shift otherwise
     warning('Number of active elements is %d, not 3. Applying zero delay to all active elements.', num_active);
     target_phases_rad = zeros(1, num_active);
else
    error('Cannot apply phases, no active elements defined.');
end
% Basic check
if length(target_phases_rad) ~= num_active
    error('Internal error: Number of phases (%d) does not match number of unique active elements (%d)', length(target_phases_rad), num_active);
end
% Calculate delays relative to the center frequency
target_delays_sec = target_phases_rad / (2 * pi * fc);
% Make delays relative (optional, but good practice)
target_delays_sec = target_delays_sec - min(target_delays_sec);
target_delays_sec_col = target_delays_sec(:); % Ensure column vector
fprintf('Calculated time delays for phase shifts.\n');

% --- Apply Time Delays using ele_delay (Corrected Indices) ---
element_no_for_handle = (1:num_active)'; % Use sequential indices [1; 2; 3; ...]
delays_col = target_delays_sec_col(:); % Ensure column vector of delays
% Check for size mismatch before calling ele_delay
if length(delays_col) ~= num_active
    error('Mismatch between number of delays (%d) and number of active elements (%d). Cannot apply ele_delay.', length(delays_col), num_active);
end
% Call ele_delay if function exists
if exist('ele_delay','file') == 2
    try
        ele_delay(Th, element_no_for_handle, delays_col);
        fprintf('Applied time delays using ele_delay to %d elements (Handle indices %s).\n', num_active, mat2str(element_no_for_handle'));
    catch ME
        warning('ele_delay failed: %s\nCheck element indices and delays. Field II internal errors often follow.', ME.message);
        fprintf('Proceeding without explicit element delays applied.\n');
    end
else
    warning('ele_delay.m not found. Cannot apply phase shifts via time delays.');
    fprintf('Proceeding without explicit element delays.\n');
end

% --- Set Impulse (Receive) and Excitation (Transmit) ---
xdc_impulse(Th, impulseResponse);    % Defines transducer receive sensitivity / transmit pulse shape component
xdc_excitation(Th, excitationPulse); % Defines the electrical excitation signal (now just '1')
fprintf('Set impulse response (Gauspuls) and excitation (Simple Delta).\n');

% --- Plot Geometry (Figure 1) ---
figure(1); clf;
% Plotting code (same as before, finding active centers)
[rows_active, cols_active] = find(enabled); active_centers_plot = zeros(num_active, 3); center_offset_x = (num_x - 1)/2 * (element_width + kerf_x); center_offset_y = (num_y - 1)/2 * (element_height + kerf_y);
for i = 1:num_active; iy = rows_active(i); ix = cols_active(i); y_pos = (iy-1)*(element_height + kerf_y) - center_offset_y; x_pos = (ix-1)*(element_width + kerf_x) - center_offset_x; active_centers_plot(i, :) = [x_pos, y_pos, 0]; end
plot3(physical_element_centers(:,1)*1e3, physical_element_centers(:,2)*1e3, physical_element_centers(:,3)*1e3, 'k.', 'MarkerSize', 1, 'DisplayName', 'Grid Elements'); hold on; plot3(active_centers_plot(:,1)*1e3, active_centers_plot(:,2)*1e3, active_centers_plot(:,3)*1e3, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Enabled Tx/Rx');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
title(sprintf('pMUT Array Geometry (%d Active Elements)', num_active)); % Generic title
axis equal; grid on; view(2); legend('Location', 'best'); set(gcf, 'Color', 'w');
saveas(gcf, 'Merged_SimplePulse_pMUT_Geometry.png');

%% Generate Measurement Matrix H using calc_hhp
fprintf('Calculating Measurement Matrix H using calc_hhp...\n');
hydrophone_positions = [X(:), Y(:), Z(:)];
if N_pixels ~= size(hydrophone_positions, 1); error('N_pixels definition mismatch'); end
[H, start_time] = calc_hhp(Th, Th, hydrophone_positions); % Simulation happens here
K = size(H, 1); M = K; N = N_pixels; % Define measurement dimensions
if N ~= size(H, 2); error('Mismatch between N_pixels (%d) and columns of H (%d)', N, size(H, 2)); end
compression = N/M;
fprintf('Measurement Matrix H (Simple Pulse/Phase): %d rows x %d columns\n', M, N); % Updated description
fprintf('Compression Ratio (N/M) = %g\n', compression);

% --- Plot Waveforms (Figure 3) ---
figure(3); clf; hold on; indices_to_plot = round(linspace(1, N, min(N, 10)));
t_axis_plot = (0:K-1)/fs * 1e6 + start_time*1e6; % Time axis in microseconds
for n = indices_to_plot; plot( t_axis_plot, H(:, n), 'DisplayName', sprintf('Pixel %d', n)); end
hold off; xlabel('Time (µs)'); ylabel('Amplitude');
title('Pulse-Echo Waveforms (Columns of H - Simple Pulse/Phase)'); % Updated title
axis tight; grid on; legend('Location', 'best'); set(gcf, 'Color', 'w');
saveas(gcf, 'Merged_SimplePulse_PulseEchoWaveforms.png');

%% Generate Scene (v) - Point targets
fprintf('Generating scene...\n');
xmin = min(x_grid); xmax = max(x_grid); zmin = min(z_grid); zmax = max(z_grid); zrange = zmax - zmin;
targets.x = [xmin*0.5, 0, xmax*0.5, 0]; targets.y = zeros(size(targets.x)); targets.z = [zmin + 0.2*zrange, zmin + 0.4*zrange, zmin + 0.6*zrange, zmin + 0.8*zrange]; amplitudes = transpose([1, 0.8, 1.2, 0.9]);
scene = zeros(length(z_grid), length(x_grid));
for i = 1:length(targets.x); [~, ix] = min(abs(x_grid - targets.x(i))); [~, iz] = min(abs(z_grid - targets.z(i))); if ix > 0 && ix <= length(x_grid) && iz > 0 && iz <= length(z_grid); scene(iz, ix) = amplitudes(i); else warning('Target %d at (%.2f, %.2f) is outside the defined grid.', i, targets.x(i)*1e3, targets.z(i)*1e3); end; end;
v = scene(:); % Ground truth image vector

% --- Plot Scene (Figure 4) ---
figure(4); clf; subplot(1, 2, 1); imagesc(x_grid*1e3, z_grid*1e3, scene); axis image; colormap gray; colorbar; xlabel('x (mm)'); ylabel('z (mm)'); title('v (matrix form)'); set(gca, 'YDir','normal'); subplot(1, 2, 2); stem(v); axis square tight; title('v (vector form)'); xlabel('Pixel Index'); ylabel('Amplitude'); sgtitle('True Image'); set(gcf, 'Color', 'w'); set(gcf, 'Position', [100 100 800 400]);
saveas(gcf, 'Merged_SimplePulse_TrueImage.png');

%% Image formation model: u = H*v + n - (Includes Debugging)
fprintf('Simulating measurement acquisition...\n');
Hv = H * v; % Calculate the ideal signal (measurements without noise)

% --- Debugging Hv ---
fprintf('Debug: Size of Hv is %s\n', mat2str(size(Hv)));
fprintf('Debug: Max absolute value in Hv = %g\n', max(abs(Hv(:))));
fprintf('Debug: Mean absolute value in Hv = %g\n', mean(abs(Hv(:))));
fprintf('Debug: Does Hv contain NaN? %d\n', any(isnan(Hv(:))));
calculated_signal_power = mean(Hv(:).^2); % Use Hv(:) to ensure vector mean
fprintf('Debug: Calculated signal_power = %g (eps = %g)\n', calculated_signal_power, eps);
% --- End Debugging ---

% Add Gaussian noise
rng('default'); % For reproducible noise
electronic_SNR_db = 30; % Define desired SNR in dB
electronic_SNR = 10^(electronic_SNR_db/10);
signal_power = calculated_signal_power; % Use the debugged value

% Check signal power before calculating noise variance
if isnan(signal_power) || signal_power < 1e-20 % Use a small threshold, maybe safer than eps
    warning('Signal power is effectively zero or NaN (%g). Cannot set noise based on SNR. Using small default noise sigma.', signal_power);
    noise_sigma = 1e-9; % Arbitrary small noise level
    actual_SNR_db = -Inf; % Indicate SNR calculation failed
else
    % Calculate noise based on signal power and desired SNR
    noise_variance = signal_power / electronic_SNR;
    noise_sigma = sqrt(noise_variance);
    actual_SNR_db = 10*log10(signal_power / noise_variance); % Actual SNR based on calculated power
end

n = noise_sigma * randn(size(Hv)); % Generate noise
u = Hv + n; % Measured signal (RF data)
fprintf('Added Gaussian noise. Target SNR: %.1f dB. Actual SNR: %.1f dB. Noise sigma: %g\n', electronic_SNR_db, actual_SNR_db, noise_sigma);

% --- Plot Measurements (Figure 6) ---
figure(6); clf; t_axis_plot = (0:K-1)/fs * 1e6 + start_time*1e6; % Time axis in microseconds
subplot(1, 3, 1); plot(t_axis_plot, Hv); axis tight; grid on; xlabel('Time (µs)'); title('Measurements (Hv)'); ylabel('Amplitude'); yl = ylim; legend('Hv');
subplot(1, 3, 2); plot(t_axis_plot, n); axis tight; grid on; xlabel('Time (µs)'); title(sprintf('Noise (SNR \\approx %.1f dB)', actual_SNR_db)); ylabel('Amplitude'); ylim(yl); legend('n'); % Show actual calculated SNR
subplot(1, 3, 3); plot(t_axis_plot, u); axis tight; grid on; xlabel('Time (µs)'); title('Measurements + Noise (u)'); ylabel('Amplitude'); ylim(yl); legend('u = Hv + n');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [100 100 1200 400]);
saveas(gcf, 'Merged_SimplePulse_Hv+n.png');

%% Reconstruction - USING ALGORITHMS FROM ORIGINAL SCRIPT (Ensemble.m)
% This section uses the H matrix, the measurement vector u (b), and noise_sigma
% derived from the simple pulse simulation above.
fprintf('Starting Reconstruction using original algorithms...\n');

% Define variables consistent with original script's reconstruction section
A = H;         % System matrix
At = transpose(A); % Transposed system matrix
b = u;         % Measurement vector (u = Hv + n)
I = scene;     % True image (for comparison)
imageResolution = size(scene); % For ADMM reshaping

% Define function handles needed by original PCG and Pinv
Afun = @(x) A*x;
Atfun = @(x) At*x;

% --- Least Norm solution with Predetermined Conjugate Gradient (PCG) ---
% --- FROM ORIGINAL SCRIPT ---
fprintf('  Calculating Least Norm (PCG) Solution (Original Method)...\n');
tic
% Define AAtfun needed specifically for this PCG formulation
AAtfun  = @(x) reshape(Afun(Atfun( x )), [M 1]); % Original way: apply A*(At*x)
maxItersCG = 1000; % Max iterations for PCG
% Use noise_sigma (calculated based on signal power if possible) as tolerance
tolCG = noise_sigma;
% Run PCG
if M > 0 && N > 0
    [x_cg, flag_cg, relres_cg, iter_cg, resvec_cg] = pcg(AAtfun, b(:), tolCG, maxItersCG);
    if flag_cg ~= 0; warning('PCG did not converge (Flag=%d). Result might be inaccurate.', flag_cg); end
    fprintf('    PCG: Flag=%d, Iters=%d, RelRes=%g\n', flag_cg, iter_cg, relres_cg);
    x_pcg = Atfun(x_cg); % Get image estimate from intermediate vector
else
    fprintf('    Skipping PCG (Zero Dim M=%d, N=%d).\n', M, N);
    x_pcg = zeros(N,1); iter_cg = 0; relres_cg = NaN; flag_cg = -1; % Placeholders
end
% Post-process and normalize result
x_pcg = real(x_pcg); x_pcg = x_pcg - min(x_pcg(:)); if max(x_pcg(:)) > 0; x_pcg = x_pcg ./ max(x_pcg(:)); end
x_pcg2D = reshape(x_pcg, imageResolution); % Reshape into 2D image
% Calculate metrics against normalized true image
norm_I = I(:); if max(norm_I) > 0; norm_I = norm_I ./ max(norm_I); end
MSE_pcg = mean( (x_pcg(:) - norm_I).^2 ); PSNR_pcg = 10*log10(1/MSE_pcg); runtime_pcg = toc;
fprintf('  Original PCG: MSE = %g, PSNR = %g dB, Time = %g s\n', MSE_pcg, PSNR_pcg, runtime_pcg);

% --- Least Norm solution with Moore-Penrose Pseudo-inverse ---
% --- FROM ORIGINAL SCRIPT ---
fprintf('  Calculating Least Norm (Pseudo-inverse) Solution (Original Method)...\n');
tic
tolPinv = 1e-6; % Tolerance for pinv
x_pinv = zeros(N, 1); % Initialize
pinv_method = 'Skipped';
if M > 0 && N > 0
    try
        AAt = A*At;
        if rank(AAt) < min(size(AAt)); fprintf('    Warning: A*At is rank deficient. Pinv may be slow or unstable.\n'); end
        AAtinv = pinv(AAt, tolPinv); % Calculate pinv(A*At)
        x_pinv = Atfun(AAtinv * b(:)); % Apply At * pinv(A*At) * b
        pinv_method = 'At*pinv(A*At)*b';
        fprintf('    Used Pinv method: %s\n', pinv_method);
    catch ME
        fprintf('    Pinv(A*At) failed: %s\n', ME.message);
        pinv_method = 'Failed';
    end
else
     fprintf('    Skipping Pinv (Zero Dim M=%d, N=%d).\n', M, N);
     pinv_method = 'Skipped (Zero Dim)';
end
% Post-process and normalize result
x_pinv = real(x_pinv); x_pinv = x_pinv - min(x_pinv(:)); if max(x_pinv(:)) > 0; x_pinv = x_pinv ./ max(x_pinv(:)); end
x_pinv2D = reshape(x_pinv, imageResolution);
% Calculate metrics
MSE_pinv = mean( (x_pinv(:) - norm_I).^2 ); PSNR_pinv = 10*log10(1/MSE_pinv); runtime_pinv = toc;
fprintf('  Original Pinv (%s): MSE = %g, PSNR = %g dB, Time = %g s\n', pinv_method, MSE_pinv, PSNR_pinv, runtime_pinv);

% --- Plot Least Norm Results (Figure 7) ---
figure(7); clf;
subplot(1, 3, 1); imagesc(x_grid*1e3, z_grid*1e3, I); axis image; colormap gray; colorbar; set(gca, 'YDir','normal'); title('Ground Truth'); xlabel('x (mm)'); ylabel('z (mm)');
subplot(1, 3, 2); imagesc(x_grid*1e3, z_grid*1e3, x_pcg2D); axis image; colormap gray; colorbar; set(gca, 'YDir','normal'); title(sprintf('Least Norm (PCG Orig.)\nIters=%d, PSNR=%.2fdB', iter_cg, PSNR_pcg)); xlabel('x (mm)'); ylabel('z (mm)');
subplot(1, 3, 3); imagesc(x_grid*1e3, z_grid*1e3, x_pinv2D); axis image; colormap gray; colorbar; set(gca, 'YDir','normal'); title(sprintf('Least Norm (Pinv Orig.)\n%s, PSNR=%.2fdB', pinv_method, PSNR_pinv)); xlabel('x (mm)'); ylabel('z (mm)');
sgtitle(sprintf('Least Norm Reconstructions (pMUT Array, Simple Pulse, SNR \\approx %.1f dB)', actual_SNR_db)); % Updated title
set(gcf, 'Color', 'w'); set(gcf, 'Position', [100 100 1200 450]);
saveas(gcf, 'Merged_SimplePulse_Reconstructed_LeastNorm.png');

%% ADMM with TV regularization - USING ALGORITHM FROM ORIGINAL SCRIPT
fprintf('  Calculating ADMM-TV Solution (Original Method)...\n');
% Normalize system matrix for ADMM stability
H_norm_factor = max(abs(H(:))); if H_norm_factor < eps; H_norm_factor = 1; end
A_admm = H ./ H_norm_factor; At_admm = transpose(A_admm);
Afun_admm = @(x) A_admm*x(:);
Atfun_admm = @(y) reshape(At_admm*y, imageResolution);
AtAfun_admm = @(x) Atfun_admm(Afun_admm(x));
% Define difference operators for TV
[Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
opDx = @(x) [Dx_sparse*x(:), Dy_sparse*x(:)];
opDtx = @(v) reshape(Dx_sparse'*v(:,1) + Dy_sparse'*v(:,2), imageResolution);
opDtDx = @(x) opDtx(opDx(x));
% ADMM parameters
numItersADMM = 25; rho = 10; lambda = 0.01;
% Scale measurement vector and noise sigma for normalized system
noise_sigma_admm = noise_sigma / H_norm_factor;
b_admm = b / H_norm_factor;
% Initialize ADMM variables
x = zeros(imageResolution); z = zeros([prod(imageResolution) 2]); u = zeros([prod(imageResolution) 2]);
% Define function handle for PCG within ADMM (original style)
Hfun = @(x_vec) reshape( AtAfun_admm(reshape(x_vec,imageResolution)) + rho.*opDtDx(reshape(x_vec,imageResolution)), [prod(imageResolution) 1]);
% Store results per iteration
PSNR_admm = zeros([numItersADMM 1]); residuals_admm = zeros([numItersADMM 1]);
% --- ADMM Loop ---
figure(8); clf; tic
for k=1:numItersADMM
    % x update
    v = z - u; bb_admm = Atfun_admm(b_admm) + rho*opDtx(v);
    maxItersCG_admm = 25; tolCG_admm = 1e-5;
    [x_vec, flag_pcg_admm, relres_pcg_admm, iter_pcg_admm] = pcg(Hfun, bb_admm(:), tolCG_admm, maxItersCG_admm, [], [], x(:));
    if flag_pcg_admm ~= 0; fprintf('    ADMM Iter %d: Inner PCG Flag=%d, Iters=%d, RelRes=%g\n', k, flag_pcg_admm, iter_pcg_admm, relres_pcg_admm); end
    x = reshape(x_vec, imageResolution);
    % Scale for display/PSNR
    x_scaled = real(x); x_min = min(x_scaled(:)); x_max = max(x_scaled(:)); if (x_max - x_min) > eps; x_scaled = (x_scaled - x_min) / (x_max - x_min); else x_scaled = zeros(imageResolution); end

    % z update (vectorial soft thresholding)
    kappa = lambda/rho; v = opDx(x) + u;
    v_norm = sqrt(sum(v.^2, 2)); v_norm(v_norm < eps) = eps;
    shrinkage_factor = max(0, 1 - kappa ./ v_norm); z = v .* shrinkage_factor;

    % u update
    u = u + opDx(x) - z;

    % Metrics
    r1 = b_admm - Afun_admm(x); r2 = opDx(x); residuals_admm(k) = 0.5*sum(r1(:).^2) + lambda.*sum( sqrt(sum(r2.^2, 2)) );
    MSE_admm = mean( (x_scaled(:) - norm_I).^2 ); PSNR_admm(k) = 10*log10(1/MSE_admm);
    runtime_ADMM = toc;

    % Plotting update
    if mod(k, 5) == 0 || k == 1 || k == numItersADMM
        subplot(2,3,[1 4]); imagesc(x_grid*1e3, z_grid*1e3, I); axis image; colormap(gca, gray); colorbar; set(gca,'YDir','normal'); title('Target Image'); xlabel('x (mm)'); ylabel('z (mm)');
        subplot(2,3,[2 5]); imagesc(x_grid*1e3, z_grid*1e3, x_scaled); axis image; colormap(gca, gray); colorbar; set(gca,'YDir','normal'); title(sprintf('ADMM TV Iter %d\nPSNR = %.2f dB', k, PSNR_admm(k))); xlabel('x (mm)'); ylabel('z (mm)');
        subplot(2,3,3); plot(1:k, PSNR_admm(1:k), 'LineWidth', 2, 'color', [1 0 1]); title('PSNR'); xlabel('Iteration'); ylabel('PSNR (dB)'); grid on; axis tight; if k > 1; ylim([min(PSNR_admm(1:k))-1, max(PSNR_admm(1:k))+1]); end
        subplot(2,3,6); plot(1:k, log10(residuals_admm(1:k)), 'LineWidth', 2); title('log10(Objective)'); xlabel('Iteration'); ylabel('Value'); grid on; axis tight; if k > 1; ylim([min(log10(residuals_admm(1:k)))-0.5, max(log10(residuals_admm(1:k)))+0.5]); end
        drawnow;
    end
end % End ADMM iterations loop
runtime_ADMM = toc;
fprintf('  Original ADMM-TV complete: Final PSNR = %g dB, Time = %gs\n', PSNR_admm(end), runtime_ADMM);
sgtitle(sprintf('ADMM TV (Original Alg., pMUT Array, Simple Pulse, SNR \\approx %.1f dB)', actual_SNR_db)); % Updated title
set(gcf, 'Color', 'w'); set(gcf, 'Position', [100 100 1200 700]);
saveas(gcf, 'Merged_SimplePulse_ADMM_TV.png');

%% End Field II Simulation
if exist('field_end','file')
    field_end;
    disp('Field II ended.');
end

%% --- Helper function for TV operators (Needed by Original ADMM Code) ---
function [Dx, Dy] = createDifferenceOperators(imageSize)
    % Creates sparse difference matrices Dx and Dy for TV regularization
    rows = imageSize(1); cols = imageSize(2); N = rows * cols;
    Dx = spdiags([-ones(N,1), ones(N,1)], [0, rows], N, N);
    % Boundary condition: Zero difference crossing columns at the right edge
    Dx(N-rows+1:N, :) = 0;
    Dy = spdiags([-ones(N,1), ones(N,1)], [0, 1], N, N);
    % Boundary condition: Zero difference crossing rows at the bottom edge
    idx = reshape(1:N, rows, cols); idx_bottom = idx(rows, :);
    Dy(idx_bottom, :) = 0;
end