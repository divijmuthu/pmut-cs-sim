%% Fixed ADMM Reconstruction with Proper TV Proximal Operator
% Fix the ADMM convergence issues with proper TV implementation

clear; clc; close all;

%% ===== CONFIGURATION =====
fprintf('=== FIXED ADMM RECONSTRUCTION ===\n');

% Create output directory
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = sprintf('fixed_admm_output/%s', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Saving results to: %s\n', output_folder);

% Test configuration
config = struct();
config.description = 'Fixed ADMM Test: 5mm targets, 20mm spacing';
config.target_size_mm = 5.0;
config.grid_spacing_mm = 20.0;
config.lambda_tv_reg = 0.1;

fprintf('Configuration: %s\n', config.description);

%% ===== STEP 1: CREATE TARGET SCENE =====
fprintf('\n=== STEP 1: CREATING TARGET SCENE ===\n');
[target_scene, target_positions, imaging_grid] = create_synthetic_target_scene(config);

%% ===== STEP 2: GENERATE SYNTHETIC H MATRIX =====
fprintf('\n=== STEP 2: GENERATING SYNTHETIC H MATRIX ===\n');
[H_matrix, transducer_positions] = generate_synthetic_h_matrix(imaging_grid);

%% ===== STEP 3: SIMULATE MEASUREMENTS =====
fprintf('\n=== STEP 3: SIMULATING MEASUREMENTS ===\n');
[measurements, measurement_times] = simulate_measurements(target_scene, H_matrix);

%% ===== STEP 4: FIXED ADMM RECONSTRUCTION =====
fprintf('\n=== STEP 4: FIXED ADMM RECONSTRUCTION ===\n');

% Test different parameter sets
test_params = [
    0.1, 1.0, 100;   % lambda_tv, rho, max_iter
    0.5, 1.0, 100;
    1.0, 1.0, 100;
    0.1, 2.0, 100;
    0.5, 2.0, 100;
    1.0, 2.0, 100;
];

results = struct();

for test_idx = 1:size(test_params, 1)
    lambda_tv = test_params(test_idx, 1);
    rho = test_params(test_idx, 2);
    max_iter = test_params(test_idx, 3);
    
    fprintf('Testing: lambda_tv=%.2f, rho=%.1f, max_iter=%d\n', lambda_tv, rho, max_iter);
    
    % Create test config
    test_config = config;
    test_config.lambda_tv_reg = lambda_tv;
    
    % Test parameters
    params = struct();
    params.rho_admm = rho;
    params.admm_max_iter = max_iter;
    params.admm_tol = 1e-6;
    
    % Perform fixed reconstruction
    [reconstructed_image, reconstruction_time, reconstruction_info] = ...
        perform_fixed_admm_reconstruction(measurements, H_matrix, test_config, params, imaging_grid);
    
    % Calculate metrics
    metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image);
    
    % Store results
    results(test_idx).lambda_tv = lambda_tv;
    results(test_idx).rho = rho;
    results(test_idx).max_iter = max_iter;
    results(test_idx).psnr = metrics.psnr;
    results(test_idx).correlation = metrics.correlation;
    results(test_idx).relative_error = metrics.relative_error;
    results(test_idx).reconstruction_time = reconstruction_time;
    results(test_idx).iterations = reconstruction_info.iterations;
    results(test_idx).converged = reconstruction_info.converged;
    results(test_idx).final_residual = reconstruction_info.final_primal_residual;
    
    fprintf('  Results: PSNR=%.2f dB, Corr=%.4f, Iter=%d, Conv=%d\n', ...
        metrics.psnr, metrics.correlation, reconstruction_info.iterations, reconstruction_info.converged);
    
    % Save individual result
    save(fullfile(output_folder, sprintf('result_%d.mat', test_idx)), ...
        'reconstructed_image', 'target_scene', 'H_matrix', 'measurements', ...
        'metrics', 'reconstruction_info', 'test_config', 'params');
end

%% ===== STEP 5: ANALYZE RESULTS =====
fprintf('\n=== STEP 5: ANALYZING RESULTS ===\n');
analyze_fixed_results(results, output_folder);

%% ===== HELPER FUNCTIONS =====

function [target_scene, target_positions, imaging_grid] = create_synthetic_target_scene(config)
    % Create synthetic target scene with high contrast
    
    % Grid parameters
    grid_width_m = 0.1;
    grid_step_m = 0.0025;
    target_distance_m = 0.05;
    grid_depth_range_m = 0.1;
    
    % Calculate imaging grid
    x_coords_img = -grid_width_m/2 : grid_step_m : grid_width_m/2;
    z_coords_img = (target_distance_m - grid_depth_range_m/2) : grid_step_m : (target_distance_m + grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    
    % Initialize target scene with background
    target_scene = ones(size(X_mesh)) * 0.1;  % Background intensity
    target_positions = [];
    
    % Create 3x3 grid of targets
    grid_center_x = 0;
    grid_center_z = target_distance_m;
    
    % Target variations for challenge - HIGH CONTRAST
    target_sizes = [config.target_size_mm * 0.8, config.target_size_mm, config.target_size_mm * 1.2];
    target_intensities = [0.9, 1.0, 1.1];  % High contrast relative to background
    
    target_idx = 1;
    for row = -1:1
        for col = -1:1
            % Calculate target position (convert to meters)
            target_x = grid_center_x + col * (config.grid_spacing_mm/1000);
            target_z = grid_center_z + row * (config.grid_spacing_mm/1000);
            
            % Select target size and intensity
            size_idx = mod(target_idx-1, length(target_sizes)) + 1;
            intensity_idx = mod(target_idx-1, length(target_intensities)) + 1;
            
            target_size = target_sizes(size_idx) / 1000;  % Convert to meters
            target_intensity = target_intensities(intensity_idx);
            
            % Create target (Gaussian profile with HIGH CONTRAST)
            target_radius_pixels = max(1, round(target_size / grid_step_m));
            target_distance = sqrt((X_mesh - target_x).^2 + (Z_mesh - target_z).^2);
            
            % Create high-contrast target (much higher than background)
            target_profile = target_intensity * exp(-(target_distance / (target_radius_pixels * 0.5)).^2);
            
            % Add target to scene (targets are much brighter than background)
            target_scene = target_scene + target_profile;
            
            % Store target position (in meters for Field II)
            target_positions = [target_positions; target_x, target_z, target_size*1000, target_intensity];
            
            target_idx = target_idx + 1;
        end
    end
    
    % Create imaging grid structure
    imaging_grid = struct('X', X_mesh, 'Z', Z_mesh, 'x_grid', x_coords_img, 'z_grid', z_coords_img, 'N_pixels', numel(X_mesh));
    
    fprintf('  Created %d targets with sizes %.1f-%.1f mm and intensities %.1f-%.1f\n', ...
        size(target_positions, 1), min(target_positions(:,3)), max(target_positions(:,3)), ...
        min(target_positions(:,4)), max(target_positions(:,4)));
    fprintf('  Target scene stats: min=%.4g, max=%.4g, mean=%.4g\n', ...
        min(target_scene(:)), max(target_scene(:)), mean(target_scene(:)));
end

function [H_matrix, transducer_positions] = generate_synthetic_h_matrix(imaging_grid)
    % Generate synthetic H matrix with realistic properties
    
    num_acquisitions = 20;
    num_pixels = imaging_grid.N_pixels;
    
    fprintf('  Generating synthetic H matrix: %d x %d\n', num_acquisitions, num_pixels);
    
    % Create synthetic H matrix with realistic properties
    % Use random matrix with controlled coherence
    rng(42); % For reproducibility
    
    % Generate base random matrix
    H_base = randn(num_acquisitions, num_pixels);
    
    % Apply spatial correlation to simulate realistic ultrasound physics
    % Each acquisition has a different spatial pattern
    for i = 1:num_acquisitions
        % Create spatial pattern based on acquisition position
        acquisition_angle = (i-1) * 2*pi / num_acquisitions;
        
        % Reshape to 2D for spatial operations
        H_2d = reshape(H_base(i, :), length(imaging_grid.x_grid), length(imaging_grid.z_grid));
        
        % Apply spatial correlation based on acquisition angle
        [X, Z] = meshgrid(imaging_grid.x_grid, imaging_grid.z_grid);
        spatial_weight = cos(acquisition_angle) * X + sin(acquisition_angle) * Z;
        H_2d = H_2d .* (1 + 0.5 * spatial_weight);
        
        % Reshape back to 1D
        H_base(i, :) = H_2d(:);
    end
    
    % Normalize columns
    column_norms = sqrt(sum(H_base.^2, 1));
    H_matrix = H_base ./ column_norms;
    
    % Add slight noise to break perfect correlations
    noise_level = 0.01;
    H_matrix = H_matrix + noise_level * randn(size(H_matrix));
    
    % Re-normalize
    column_norms = sqrt(sum(H_matrix.^2, 1));
    H_matrix = H_matrix ./ column_norms;
    
    % Calculate transducer positions (synthetic)
    transducer_positions = linspace(-0.04, 0.04, 10);
    
    fprintf('  H matrix generation completed\n');
    fprintf('  - Matrix dimensions: %d x %d\n', size(H_matrix, 1), size(H_matrix, 2));
end

function [measurements, measurement_times] = simulate_measurements(target_scene, H_matrix)
    % Simulate measurements
    
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

function [reconstructed_image, reconstruction_time, reconstruction_info] = perform_fixed_admm_reconstruction(measurements, H_matrix, config, params, imaging_grid)
    % Perform ADMM reconstruction with FIXED TV proximal operator
    
    tic;
    
    % Calculate imaging grid dimensions
    num_x = length(imaging_grid.x_grid);
    num_z = length(imaging_grid.z_grid);
    
    % Initialize reconstruction
    num_pixels = size(H_matrix, 2);
    
    % ADMM parameters
    lambda_tv = config.lambda_tv_reg;
    rho = params.rho_admm;
    max_iter = params.admm_max_iter;
    tol = params.admm_tol;
    
    fprintf('    ADMM parameters: lambda_tv=%.2f, rho=%.1f, max_iter=%d, tol=%.1e\n', ...
        lambda_tv, rho, max_iter, tol);
    
    % Initialize variables
    x = zeros(num_pixels, 1);  % Main variable
    z = zeros(num_pixels, 1);  % Auxiliary variable
    u = zeros(num_pixels, 1);  % Dual variable
    
    % Precompute matrices
    HtH = H_matrix' * H_matrix;
    Hty = H_matrix' * measurements;
    
    % ADMM iterations with FIXED TV proximal operator
    primal_residuals = zeros(max_iter, 1);
    dual_residuals = zeros(max_iter, 1);
    objective_values = zeros(max_iter, 1);
    
    fprintf('    Starting FIXED ADMM iterations...\n');
    
    for iter = 1:max_iter
        % x-update (least squares with regularization)
        x_prev = x;
        x = (HtH + rho * eye(num_pixels)) \ (Hty + rho * (z - u));
        
        % z-update (FIXED TV proximal operator)
        z_prev = z;
        z = fixed_prox_tv(x + u, lambda_tv / rho, num_x, num_z);
        
        % u-update (dual variable)
        u = u + (x - z);
        
        % Calculate residuals and objective
        primal_residual = norm(x - z);
        dual_residual = rho * norm(z - z_prev);
        
        primal_residuals(iter) = primal_residual;
        dual_residuals(iter) = dual_residual;
        
        % Calculate objective value
        data_term = 0.5 * norm(H_matrix * x - measurements)^2;
        tv_term = lambda_tv * fixed_tv_norm(x, num_x, num_z);
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
    
    % Store reconstruction info
    reconstruction_info.iterations = iter;
    reconstruction_info.converged = (primal_residual < tol && dual_residual < tol);
    reconstruction_info.final_primal_residual = primal_residual;
    reconstruction_info.final_dual_residual = dual_residual;
    reconstruction_info.final_objective = objective_values(iter);
    reconstruction_info.primal_residuals = primal_residuals(1:iter);
    reconstruction_info.dual_residuals = dual_residuals(1:iter);
    reconstruction_info.objective_values = objective_values(1:iter);
    
    % Reshape to image
    reconstructed_image = reshape(x, num_x, num_z);
    
    fprintf('    FIXED ADMM completed in %.2f seconds\n', reconstruction_time);
    fprintf('    - Iterations: %d\n', iter);
    fprintf('    - Converged: %s\n', mat2str(reconstruction_info.converged));
    fprintf('    - Final primal residual: %.2e\n', primal_residual);
end

function z = fixed_prox_tv(x, lambda, num_x, num_z)
    % FIXED TV proximal operator using proper gradient descent
    
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

function tv_val = fixed_tv_norm(x, num_x, num_z)
    % Calculate TV norm properly
    x_reshaped = reshape(x, num_x, num_z);
    
    % Calculate gradients
    grad_x = diff(x_reshaped, 1, 1);
    grad_z = diff(x_reshaped, 1, 2);
    
    % TV norm
    tv_val = sum(abs(grad_x(:))) + sum(abs(grad_z(:)));
end

function metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image)
    % Calculate reconstruction quality metrics with meaningful comparison
    
    % Use more meaningful normalization - scale to same range but preserve structure
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
    
    % Calculate PSNR using the normalized images
    mse = mean((target_norm(:) - recon_norm(:)).^2);
    
    if mse > 0
        metrics.psnr = 10 * log10(1 / mse);
    else
        metrics.psnr = Inf;
    end
    
    % Calculate correlation
    correlation_matrix = corrcoef(target_norm(:), recon_norm(:));
    metrics.correlation = correlation_matrix(1, 2);
    
    % Calculate structural similarity (if Image Processing Toolbox available)
    try
        metrics.ssim = ssim(recon_norm, target_norm);
    catch
        metrics.ssim = NaN;
    end
    
    % Calculate relative error
    metrics.relative_error = norm(target_norm(:) - recon_norm(:)) / norm(target_norm(:));
end

function analyze_fixed_results(results, output_folder)
    % Analyze fixed ADMM results
    
    fprintf('  === FIXED ADMM RESULTS ANALYSIS ===\n');
    
    % Find best results
    psnr_values = [results.psnr];
    correlation_values = [results.correlation];
    converged_values = [results.converged];
    
    [best_psnr, best_psnr_idx] = max(psnr_values);
    [best_correlation, best_correlation_idx] = max(correlation_values);
    converged_count = sum(converged_values);
    
    fprintf('  Best PSNR: %.2f dB (Test %d)\n', best_psnr, best_psnr_idx);
    fprintf('  Best Correlation: %.4f (Test %d)\n', best_correlation, best_correlation_idx);
    fprintf('  Convergence rate: %.1f%% (%d/%d converged)\n', ...
        converged_count/length(results)*100, converged_count, length(results));
    
    % Print all results
    fprintf('\n  All Results:\n');
    for i = 1:length(results)
        fprintf('  Test %d: lambda_tv=%.2f, rho=%.1f, PSNR=%.2f dB, Corr=%.4f, Conv=%d\n', ...
            i, results(i).lambda_tv, results(i).rho, results(i).psnr, results(i).correlation, results(i).converged);
    end
    
    fprintf('  === END FIXED ADMM RESULTS ANALYSIS ===\n');
end

fprintf('\n=== FIXED ADMM RECONSTRUCTION COMPLETE ===\n');
fprintf('All results saved to: %s\n', output_folder); 