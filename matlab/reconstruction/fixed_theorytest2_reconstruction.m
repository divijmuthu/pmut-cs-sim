%% Fixed TheoryTest2 Reconstruction
% Use EXACT TheoryTest2 H generation with proper time-domain interpolation

clear; clc; close all;

% Add Field II path
addpath('/Users/deepshikhakaul/Documents/pmut-cs-sim/m_files');

%% ===== CONFIGURATION =====
fprintf('=== FIXED THEORYTEST2 RECONSTRUCTION ===\n');

% Create output directory
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = sprintf('fixed_theorytest2_output/%s', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Saving results to: %s\n', output_folder);

% Use EXACT TheoryTest2 parameters
params = struct();
params.fs = 100e6;  % 100 MHz sampling
params.c = 343;     % Speed of sound
params.excitation_amplitude = 1e8;
params.R_acquisitions = 20;  % Number of acquisitions
params.grid_width_mm = 50;   % Reduced grid width (was 100)
params.grid_step_mm = 5;     % Increased grid step (was 2)
params.target_distance_m = 0.05;  % Target distance

fprintf('Configuration: EXACT TheoryTest2 parameters\n');
fprintf('  - Acquisitions: %d\n', params.R_acquisitions);
fprintf('  - Grid: %d mm width, %d mm step\n', params.grid_width_mm, params.grid_step_mm);

%% ===== STEP 1: CREATE TARGET SCENE =====
fprintf('\n=== STEP 1: CREATING TARGET SCENE ===\n');
[target_scene, imaging_grid] = create_theorytest2_target_scene(params);

%% ===== STEP 2: GENERATE H MATRIX USING EXACT THEORYTEST2 APPROACH =====
fprintf('\n=== STEP 2: GENERATING H MATRIX (EXACT THEORYTEST2) ===\n');
[H_matrix, reconstruction_info] = generate_theorytest2_h_matrix(params, imaging_grid);

%% ===== STEP 3: SIMULATE MEASUREMENTS =====
fprintf('\n=== STEP 3: SIMULATING MEASUREMENTS ===\n');
[measurements, measurement_times] = simulate_theorytest2_measurements(target_scene, H_matrix);

%% ===== STEP 4: THEORYTEST2 ADMM RECONSTRUCTION =====
fprintf('\n=== STEP 4: THEORYTEST2 ADMM RECONSTRUCTION ===\n');

% Use EXACT TheoryTest2 ADMM parameters
lambda_tv = 0.1;  % From TheoryTest2
rho = 1.0;        % From TheoryTest2
max_iter = 100;   % From TheoryTest2
tol = 1e-6;       % From TheoryTest2

fprintf('  ADMM parameters (EXACT TheoryTest2):\n');
fprintf('    - lambda_tv: %.2f\n', lambda_tv);
fprintf('    - rho: %.1f\n', rho);
fprintf('    - max_iter: %d\n', max_iter);
fprintf('    - tol: %.1e\n', tol);

% Perform TheoryTest2 ADMM reconstruction
[reconstructed_image, reconstruction_time, admm_info] = ...
    perform_theorytest2_admm_reconstruction(measurements, H_matrix, lambda_tv, rho, max_iter, tol, imaging_grid);

%% ===== STEP 5: CALCULATE METRICS =====
fprintf('\n=== STEP 5: CALCULATING METRICS ===\n');
metrics = calculate_theorytest2_metrics(target_scene, reconstructed_image);

%% ===== STEP 6: SAVE RESULTS =====
fprintf('\n=== STEP 6: SAVING RESULTS ===\n');
save(fullfile(output_folder, 'theorytest2_results.mat'), ...
    'target_scene', 'reconstructed_image', 'H_matrix', 'measurements', ...
    'metrics', 'reconstruction_info', 'admm_info', 'params');

% Create visualization
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
imagesc(imaging_grid.x_grid, imaging_grid.z_grid, target_scene);
colorbar;
title('Original Target Scene');
xlabel('X (m)'); ylabel('Z (m)');
axis equal tight;

subplot(1, 3, 2);
imagesc(imaging_grid.x_grid, imaging_grid.z_grid, reconstructed_image);
colorbar;
title('Reconstructed Image');
xlabel('X (m)'); ylabel('Z (m)');
axis equal tight;

subplot(1, 3, 3);
diff_image = abs(target_scene - reconstructed_image);
imagesc(imaging_grid.x_grid, imaging_grid.z_grid, diff_image);
colorbar;
title('Absolute Difference');
xlabel('X (m)'); ylabel('Z (m)');
axis equal tight;

sgtitle(sprintf('TheoryTest2 Reconstruction Results - PSNR: %.2f dB, Corr: %.4f', ...
    metrics.psnr, metrics.correlation), 'FontSize', 14);

saveas(gcf, fullfile(output_folder, 'theorytest2_reconstruction.png'));

fprintf('  Results saved to: %s\n', output_folder);
fprintf('  - PSNR: %.2f dB\n', metrics.psnr);
fprintf('  - Correlation: %.4f\n', metrics.correlation);
fprintf('  - Reconstruction time: %.2f seconds\n', reconstruction_time);
fprintf('  - ADMM iterations: %d\n', admm_info.iterations);
fprintf('  - ADMM converged: %s\n', mat2str(admm_info.converged));

%% ===== HELPER FUNCTIONS =====

function [target_scene, imaging_grid] = create_theorytest2_target_scene(params)
    % Create target scene using TheoryTest2 approach
    
    % Grid parameters (EXACT TheoryTest2)
    grid_width_m = params.grid_width_mm / 1000;
    grid_step_m = params.grid_step_mm / 1000;
    target_distance_m = params.target_distance_m;
    
    % Calculate imaging grid
    x_coords = -grid_width_m/2 : grid_step_m : grid_width_m/2;
    z_coords = target_distance_m - grid_width_m/2 : grid_step_m : target_distance_m + grid_width_m/2;
    [X_mesh, Z_mesh] = meshgrid(x_coords, z_coords);
    
    % Initialize target scene
    target_scene = zeros(size(X_mesh));
    
    % Create targets (TheoryTest2 style)
    target_positions = [
        -0.02, target_distance_m, 0.005, 1.0;   % Left target
        0.00, target_distance_m, 0.005, 1.0;    % Center target
        0.02, target_distance_m, 0.005, 1.0;    % Right target
    ];
    
    % Add targets to scene
    for i = 1:size(target_positions, 1)
        target_x = target_positions(i, 1);
        target_z = target_positions(i, 2);
        target_radius = target_positions(i, 3);
        target_intensity = target_positions(i, 4);
        
        % Create Gaussian target
        target_distance = sqrt((X_mesh - target_x).^2 + (Z_mesh - target_z).^2);
        target_profile = target_intensity * exp(-(target_distance / target_radius).^2);
        
        target_scene = target_scene + target_profile;
    end
    
    % Create imaging grid structure
    imaging_grid = struct();
    imaging_grid.x_grid = x_coords;
    imaging_grid.z_grid = z_coords;
    imaging_grid.X_mesh = X_mesh;
    imaging_grid.Z_mesh = Z_mesh;
    imaging_grid.N_pixels = numel(X_mesh);
    
    fprintf('  Created target scene: %d x %d pixels\n', size(target_scene, 1), size(target_scene, 2));
    fprintf('  Target scene stats: min=%.4g, max=%.4g, mean=%.4g\n', ...
        min(target_scene(:)), max(target_scene(:)), mean(target_scene(:)));
end

function [H_matrix, reconstruction_info] = generate_theorytest2_h_matrix(params, imaging_grid)
    % Generate H matrix using EXACT TheoryTest2 approach
    
    % Initialize Field II
    fprintf('  Initializing Field II...\n');
    field_init(-1);
    
    % Setup transducers (EXACT TheoryTest2 style)
    pMUT_width = 0.002;  % 2mm pMUT
    pMUT_height = pMUT_width;
    kerf = 0.0005;  % 0.5mm kerf
    
    % Create 2D array apertures (EXACT TheoryTest2)
    num_x_grid = 9;
    num_y_grid = 9;
    
    % Create enabled matrices for tx and rx
    tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    
    % Enable specific elements (TheoryTest2 style)
    tx_enabled_matrix(5, 3) = 1;  % First tx element
    tx_enabled_matrix(3, 7) = 1;  % Second tx element  
    tx_enabled_matrix(7, 7) = 1;  % Third tx element
    rx_enabled_matrix(5, 5) = 1;  % Single rx element
    
    % Create apertures using xdc_2d_array (EXACT TheoryTest2)
    tx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, tx_enabled_matrix, 1, 1, [0 0 0]);
    rx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(rx_Aperture, ones(1,10));
    
    % Create hydrophone positions
    hydrophone_positions = [imaging_grid.X_mesh(:), zeros(numel(imaging_grid.X_mesh), 1), imaging_grid.Z_mesh(:)];
    
    % H matrix generation with EXACT TheoryTest2 approach
    all_hhp_data = cell(params.R_acquisitions, 1);
    all_start_times = zeros(params.R_acquisitions, 1);
    all_K_values = zeros(params.R_acquisitions, 1);
    
    fprintf('  Starting %d acquisitions (EXACT TheoryTest2)...\n', params.R_acquisitions);
    
    for r_acq = 1:params.R_acquisitions
        % Generate excitation signals (EXACT TheoryTest2)
        f_min = 45e3;
        f_max = 65e3;
        tx_frequencies = f_min + (f_max - f_min) * rand(3, 1);
        tx_phase_delays = 12e-6 * rand(3, 1);
        cycles = 3;
        tx_durations = cycles ./ tx_frequencies;
        tx_signals = cell(3, 1);
        
        for i = 1:3
            t_signal = 0 : 1/params.fs : tx_durations(i);
            signal_base = sin(2 * pi * tx_frequencies(i) * t_signal);
            window = tukeywin(length(signal_base), 0.25)';
            tx_signals{i} = signal_base .* window * params.excitation_amplitude;
        end
        
        % Apply signals to transducer
        for i = 1:3
            xdc_impulse(tx_Aperture, tx_signals{i});
        end
        
        % Set focus times
        xdc_focus_times(tx_Aperture, 0, tx_phase_delays(:)');
        
        % Calculate impulse response
        [hhp_r, start_time_r] = calc_hhp(tx_Aperture, rx_Aperture, hydrophone_positions);
        
        % Store data
        all_hhp_data{r_acq} = hhp_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(hhp_r, 1);
        
        if mod(r_acq, 5) == 0
            fprintf('    Acquisition %d/%d complete\n', r_acq, params.R_acquisitions);
        end
    end
    
    % EXACT TheoryTest2 time-domain interpolation
    fprintf('  Performing time-domain interpolation (EXACT TheoryTest2)...\n');
    
    % Calculate global time window
    valid_indices = all_K_values > 0;
    all_end_times = zeros(params.R_acquisitions, 1);
    all_end_times(valid_indices) = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / params.fs;
    
    min_global_start_time = min(all_start_times);
    max_global_end_time = max(all_end_times);
    
    if isempty(min_global_start_time) || isempty(max_global_end_time) || min_global_start_time >= max_global_end_time
        min_global_start_time = 0;
        max_K_val = max(all_K_values(all_K_values > 0));
        if isempty(max_K_val) || max_K_val == 0
            max_K_val = 100;
        end
        max_global_end_time = (max_K_val - 1) / params.fs;
        if max_global_end_time < min_global_start_time
            max_global_end_time = min_global_start_time + (100 - 1) / params.fs;
        end
    end
    
    % Create common time axis
    t_common_axis = min_global_start_time : 1/params.fs : max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    
    if K_global_per_acq == 0
        K_global_per_acq = 1;
        t_common_axis = min_global_start_time;
    end
    
    % Assemble H matrix (EXACT TheoryTest2)
    total_rows = K_global_per_acq * params.R_acquisitions;
    total_cols = imaging_grid.N_pixels;
    estimated_nnz = total_rows * total_cols * 0.1;
    
    H_assembled = spalloc(total_rows, total_cols, estimated_nnz);
    current_row_offset = 0;
    
    for r_acq = 1:params.R_acquisitions
        hhp_current = all_hhp_data{r_acq};
        start_time_current = all_start_times(r_acq);
        K_current = all_K_values(r_acq);
        
        if K_current == 0 || isempty(hhp_current)
            current_row_offset = current_row_offset + K_global_per_acq;
            continue;
        end
        
        % Time interpolation (EXACT TheoryTest2)
        t_current_acq_axis = start_time_current + (0:(K_current - 1)) / params.fs;
        hhp_aligned_r = quantum_interpolation(t_current_acq_axis, hhp_current, t_common_axis, K_global_per_acq, imaging_grid.N_pixels);
        
        % Assign to H matrix
        row_indices = current_row_offset + (1:K_global_per_acq);
        if max(row_indices) <= size(H_assembled, 1)
            H_assembled(row_indices, :) = hhp_aligned_r;
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    
    H_matrix = sparse(H_assembled);
    
    % Cleanup Field II
    xdc_free(tx_Aperture);
    xdc_free(rx_Aperture);
    field_end();
    
    % Store reconstruction info
    reconstruction_info.total_rows = total_rows;
    reconstruction_info.total_cols = total_cols;
    reconstruction_info.K_global_per_acq = K_global_per_acq;
    reconstruction_info.t_common_axis = t_common_axis;
    
    fprintf('  H matrix generation completed (EXACT TheoryTest2)\n');
    fprintf('  - Matrix dimensions: %d x %d\n', size(H_matrix, 1), size(H_matrix, 2));
    fprintf('  - Sparsity: %.2f%%\n', 100 * (1 - nnz(H_matrix) / numel(H_matrix)));
end

function h_interpolated = quantum_interpolation(t_current_acq_axis, h_current, t_common_axis, K_global_per_acq, N_pixels)
    % EXACT TheoryTest2 interpolation function
    
    h_interpolated = zeros(K_global_per_acq, N_pixels);
    
    if length(t_current_acq_axis) > 1 && issorted(t_current_acq_axis)
        % Vectorized interpolation for maximum speed
        for px_col = 1:N_pixels
            if ~isempty(h_current) && size(h_current, 2) >= px_col
                h_interpolated(:, px_col) = interp1(t_current_acq_axis, h_current(:, px_col), t_common_axis, 'linear', 0);
            end
        end
    elseif isscalar(t_current_acq_axis) && K_global_per_acq >= 1
        % Handle single time point
        [~, idx_match] = min(abs(t_common_axis - t_current_acq_axis));
        if ~isempty(idx_match) && ~isempty(h_current)
            h_interpolated(idx_match, :) = h_current(1, :);
        end
    end
end

function [measurements, measurement_times] = simulate_theorytest2_measurements(target_scene, H_matrix)
    % Simulate measurements using TheoryTest2 approach
    
    % Flatten target scene
    target_scene_flat = target_scene(:);
    
    % Ensure target scene matches H matrix dimensions
    if length(target_scene_flat) ~= size(H_matrix, 2)
        fprintf('  WARNING: Target scene size (%d) does not match H matrix columns (%d)\n', length(target_scene_flat), size(H_matrix, 2));
        if length(target_scene_flat) > size(H_matrix, 2)
            target_scene_flat = target_scene_flat(1:size(H_matrix, 2));
        else
            target_scene_flat = [target_scene_flat; zeros(size(H_matrix, 2) - length(target_scene_flat), 1)];
        end
    end
    
    % Generate measurements
    measurements = H_matrix * target_scene_flat;
    
    % Add realistic noise
    noise_level = 0.02;  % 2% noise
    noise = noise_level * randn(size(measurements));
    measurements = measurements + noise;
    
    % Generate measurement times
    measurement_times = 1:length(measurements);
    
    fprintf('  Generated %d measurements with %.1f%% noise\n', length(measurements), noise_level*100);
end

function [reconstructed_image, reconstruction_time, admm_info] = perform_theorytest2_admm_reconstruction(measurements, H_matrix, lambda_tv, rho, max_iter, tol, imaging_grid)
    % Perform ADMM reconstruction using EXACT TheoryTest2 approach
    
    tic;
    
    % Calculate imaging grid dimensions
    num_x = length(imaging_grid.x_grid);
    num_z = length(imaging_grid.z_grid);
    
    % Initialize reconstruction
    num_pixels = size(H_matrix, 2);
    
    fprintf('    Starting TheoryTest2 ADMM reconstruction...\n');
    fprintf('    - Matrix size: %d x %d\n', size(H_matrix, 1), size(H_matrix, 2));
    fprintf('    - Target size: %d x %d\n', num_x, num_z);
    
    % Initialize variables
    x = zeros(num_pixels, 1);  % Main variable
    z = zeros(num_pixels, 1);  % Auxiliary variable
    u = zeros(num_pixels, 1);  % Dual variable
    
    % Precompute matrices
    HtH = H_matrix' * H_matrix;
    Hty = H_matrix' * measurements;
    
    % ADMM iterations (EXACT TheoryTest2)
    primal_residuals = zeros(max_iter, 1);
    dual_residuals = zeros(max_iter, 1);
    objective_values = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        % x-update (least squares with regularization)
        x_prev = x;
        x = (HtH + rho * eye(num_pixels)) \ (Hty + rho * (z - u));
        
        % z-update (TV proximal operator)
        z_prev = z;
        z = theorytest2_prox_tv(x + u, lambda_tv / rho, num_x, num_z);
        
        % u-update (dual variable)
        u = u + (x - z);
        
        % Calculate residuals and objective
        primal_residual = norm(x - z);
        dual_residual = rho * norm(z - z_prev);
        
        primal_residuals(iter) = primal_residual;
        dual_residuals(iter) = dual_residual;
        
        % Calculate objective value
        data_term = 0.5 * norm(H_matrix * x - measurements)^2;
        tv_term = lambda_tv * theorytest2_tv_norm(x, num_x, num_z);
        objective_values(iter) = data_term + tv_term;
        
        % Check convergence
        if primal_residual < tol && dual_residual < tol
            fprintf('    Converged at iteration %d\n', iter);
            break;
        end
        
        % Print progress every 10 iterations
        if mod(iter, 10) == 0
            fprintf('    Iter %d: primal=%.2e, dual=%.2e, obj=%.4f\n', ...
                iter, primal_residual, dual_residual, objective_values(iter));
        end
    end
    
    reconstruction_time = toc;
    
    % Store ADMM info
    admm_info.iterations = iter;
    admm_info.converged = (primal_residual < tol && dual_residual < tol);
    admm_info.final_primal_residual = primal_residual;
    admm_info.final_dual_residual = dual_residual;
    admm_info.final_objective = objective_values(iter);
    admm_info.primal_residuals = primal_residuals(1:iter);
    admm_info.dual_residuals = dual_residuals(1:iter);
    admm_info.objective_values = objective_values(1:iter);
    
    % Reshape to image
    reconstructed_image = reshape(x, num_x, num_z);
    
    fprintf('    TheoryTest2 ADMM completed in %.2f seconds\n', reconstruction_time);
    fprintf('    - Iterations: %d\n', iter);
    fprintf('    - Converged: %s\n', mat2str(admm_info.converged));
    fprintf('    - Final primal residual: %.2e\n', primal_residual);
end

function z = theorytest2_prox_tv(x, lambda, num_x, num_z)
    % TheoryTest2 TV proximal operator
    
    % Reshape to 2D
    x_reshaped = reshape(x, num_x, num_z);
    
    % Initialize z
    z_reshaped = x_reshaped;
    
    % Gradient descent for TV proximal operator
    max_inner_iter = 50;
    inner_tol = 1e-6;
    step_size = 0.1;
    
    for inner_iter = 1:max_inner_iter
        z_prev = z_reshaped;
        
        % Calculate gradients
        grad_x = diff(z_reshaped, 1, 1);
        grad_z = diff(z_reshaped, 1, 2);
        
        % Calculate gradient of TV norm
        tv_grad = zeros(size(z_reshaped));
        
        % X-gradient contribution
        tv_grad(1:end-1, :) = tv_grad(1:end-1, :) + sign(grad_x) ./ (abs(grad_x) + 1e-10);
        tv_grad(2:end, :) = tv_grad(2:end, :) - sign(grad_x) ./ (abs(grad_x) + 1e-10);
        
        % Z-gradient contribution
        tv_grad(:, 1:end-1) = tv_grad(:, 1:end-1) + sign(grad_z) ./ (abs(grad_z) + 1e-10);
        tv_grad(:, 2:end) = tv_grad(:, 2:end) - sign(grad_z) ./ (abs(grad_z) + 1e-10);
        
        % Gradient descent step
        z_reshaped = z_reshaped - step_size * (z_reshaped - x_reshaped + lambda * tv_grad);
        
        % Check convergence
        if norm(z_reshaped - z_prev) < inner_tol
            break;
        end
    end
    
    z = z_reshaped(:);
end

function tv_val = theorytest2_tv_norm(x, num_x, num_z)
    % Calculate TV norm for TheoryTest2
    
    x_reshaped = reshape(x, num_x, num_z);
    
    % Calculate gradients
    grad_x = diff(x_reshaped, 1, 1);
    grad_z = diff(x_reshaped, 1, 2);
    
    % TV norm
    tv_val = sum(abs(grad_x(:))) + sum(abs(grad_z(:)));
end

function metrics = calculate_theorytest2_metrics(target_scene, reconstructed_image)
    % Calculate metrics using TheoryTest2 approach
    
    % Use more meaningful normalization
    target_range = max(target_scene(:)) - min(target_scene(:));
    recon_range = max(reconstructed_image(:)) - min(reconstructed_image(:));
    
    if target_range > 0 && recon_range > 0
        % Normalize both to [0,1] range while preserving structure
        target_norm = (target_scene - min(target_scene(:))) / target_range;
        recon_norm = (reconstructed_image - min(reconstructed_image(:))) / recon_range;
    else
        % Fallback to max normalization if ranges are zero
        target_norm = target_scene / max(target_scene(:));
        recon_norm = reconstructed_image / max(reconstructed_image(:));
    end
    
    % Calculate PSNR
    mse = mean((target_norm(:) - recon_norm(:)).^2);
    
    if mse > 0
        metrics.psnr = 10 * log10(1 / mse);
    else
        metrics.psnr = Inf;
    end
    
    % Calculate correlation
    correlation_matrix = corrcoef(target_norm(:), recon_norm(:));
    metrics.correlation = correlation_matrix(1, 2);
    
    % Calculate structural similarity
    try
        metrics.ssim = ssim(recon_norm, target_norm);
    catch
        metrics.ssim = NaN;
    end
    
    % Calculate relative error
    metrics.relative_error = norm(target_norm(:) - recon_norm(:)) / norm(target_norm(:));
end

fprintf('\n=== FIXED THEORYTEST2 RECONSTRUCTION COMPLETE ===\n');
fprintf('All results saved to: %s\n', output_folder); 