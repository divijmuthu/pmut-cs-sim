%% ADMM Diagnostics with Synthetic H Matrix
% Test ADMM parameters using synthetic data to avoid Field II crashes

clear; clc; close all;

%% ===== CONFIGURATION =====
fprintf('=== ADMM DIAGNOSTICS WITH SYNTHETIC H MATRIX ===\n');

% Create output directory
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = sprintf('admm_diagnostics_output/%s', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Saving results to: %s\n', output_folder);

% Test configuration
config = struct();
config.description = 'Synthetic Test: 5mm targets, 20mm spacing';
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

%% ===== STEP 3: H MATRIX DIAGNOSTICS =====
fprintf('\n=== STEP 3: H MATRIX DIAGNOSTICS ===\n');
perform_h_matrix_diagnostics(H_matrix, target_scene, output_folder);

%% ===== STEP 4: SIMULATE MEASUREMENTS =====
fprintf('\n=== STEP 4: SIMULATING MEASUREMENTS ===\n');
[measurements, measurement_times] = simulate_measurements(target_scene, H_matrix);

%% ===== STEP 5: ADMM PARAMETER SWEEP =====
fprintf('\n=== STEP 5: ADMM PARAMETER SWEEP ===\n');

% Define parameter ranges to test
lambda_tv_values = [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0];
rho_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
max_iter_values = [50, 100, 200, 500];

% Initialize results storage
results = struct();
result_idx = 1;

% Test all combinations
for lambda_idx = 1:length(lambda_tv_values)
    for rho_idx = 1:length(rho_values)
        for iter_idx = 1:length(max_iter_values)
            
            lambda_tv = lambda_tv_values(lambda_idx);
            rho = rho_values(rho_idx);
            max_iter = max_iter_values(iter_idx);
            
            fprintf('Testing: lambda_tv=%.2f, rho=%.1f, max_iter=%d\n', lambda_tv, rho, max_iter);
            
            % Create test config
            test_config = config;
            test_config.lambda_tv_reg = lambda_tv;
            
            % Test parameters
            test_params = struct();
            test_params.rho_admm = rho;
            test_params.admm_max_iter = max_iter;
            test_params.admm_tol = 1e-6;
            
            % Perform reconstruction
            [reconstructed_image, reconstruction_time, reconstruction_info] = ...
                perform_admm_reconstruction_with_diagnostics(measurements, H_matrix, test_config, test_params, imaging_grid);
            
            % Calculate metrics
            metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image);
            
            % Store results
            results(result_idx).lambda_tv = lambda_tv;
            results(result_idx).rho = rho;
            results(result_idx).max_iter = max_iter;
            results(result_idx).psnr = metrics.psnr;
            results(result_idx).correlation = metrics.correlation;
            results(result_idx).relative_error = metrics.relative_error;
            results(result_idx).reconstruction_time = reconstruction_time;
            results(result_idx).iterations = reconstruction_info.iterations;
            results(result_idx).converged = reconstruction_info.converged;
            results(result_idx).final_residual = reconstruction_info.final_primal_residual;
            
            fprintf('  Results: PSNR=%.2f dB, Corr=%.4f, Iter=%d, Conv=%d\n', ...
                metrics.psnr, metrics.correlation, reconstruction_info.iterations, reconstruction_info.converged);
            
            result_idx = result_idx + 1;
        end
    end
end

%% ===== STEP 6: ANALYZE RESULTS =====
fprintf('\n=== STEP 6: ANALYZING RESULTS ===\n');
analyze_admm_results(results, output_folder);

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

function perform_h_matrix_diagnostics(H_matrix, target_scene, output_folder)
    % Perform comprehensive H matrix diagnostics
    
    fprintf('  === H MATRIX DIAGNOSTICS ===\n');
    
    % Basic matrix properties
    [m, n] = size(H_matrix);
    fprintf('  Matrix dimensions: %d x %d\n', m, n);
    fprintf('  Condition number: %.2e\n', cond(H_matrix));
    fprintf('  Rank: %d\n', rank(H_matrix));
    
    % Singular value analysis
    [U, S, V] = svd(H_matrix);
    singular_values = diag(S);
    fprintf('  Singular values: min=%.2e, max=%.2e, ratio=%.2e\n', ...
        min(singular_values), max(singular_values), max(singular_values)/min(singular_values));
    
    % Coherence analysis
    coherence_matrix = abs(H_matrix' * H_matrix);
    coherence_matrix = coherence_matrix - diag(diag(coherence_matrix)); % Remove diagonal
    max_coherence = max(coherence_matrix(:));
    mean_coherence = mean(coherence_matrix(:));
    fprintf('  Coherence: max=%.4f, mean=%.4f\n', max_coherence, mean_coherence);
    
    % Sparsity analysis
    sparsity = 1 - nnz(H_matrix) / numel(H_matrix);
    fprintf('  Sparsity: %.2f%%\n', sparsity * 100);
    
    % Column norms analysis
    column_norms = sqrt(sum(H_matrix.^2, 1));
    fprintf('  Column norms: min=%.4f, max=%.4f, mean=%.4f\n', ...
        min(column_norms), max(column_norms), mean(column_norms));
    
    % Row norms analysis
    row_norms = sqrt(sum(H_matrix.^2, 2));
    fprintf('  Row norms: min=%.4f, max=%.4f, mean=%.4f\n', ...
        min(row_norms), max(row_norms), mean(row_norms));
    
    % Test forward model with target scene
    target_flat = target_scene(:);
    if length(target_flat) ~= n
        fprintf('  WARNING: Target scene size (%d) does not match H matrix columns (%d)\n', length(target_flat), n);
        target_flat = target_flat(1:n); % Truncate if needed
    end
    
    simulated_measurements = H_matrix * target_flat;
    fprintf('  Simulated measurements: min=%.4f, max=%.4f, mean=%.4f\n', ...
        min(simulated_measurements), max(simulated_measurements), mean(simulated_measurements));
    
    % Create diagnostic plots
    figure('Position', [100, 100, 1200, 800]);
    
    % H matrix visualization
    subplot(2, 3, 1);
    imagesc(H_matrix);
    colorbar;
    title('H Matrix');
    xlabel('Pixel Index'); ylabel('Acquisition Index');
    
    % Singular values
    subplot(2, 3, 2);
    semilogy(singular_values, 'b-', 'LineWidth', 2);
    title('Singular Values');
    xlabel('Index'); ylabel('Value');
    grid on;
    
    % Coherence matrix
    subplot(2, 3, 3);
    imagesc(coherence_matrix);
    colorbar;
    title('Coherence Matrix (off-diagonal)');
    xlabel('Pixel Index'); ylabel('Pixel Index');
    
    % Column norms
    subplot(2, 3, 4);
    plot(column_norms, 'b-', 'LineWidth', 2);
    title('Column Norms');
    xlabel('Pixel Index'); ylabel('Norm');
    grid on;
    
    % Row norms
    subplot(2, 3, 5);
    plot(row_norms, 'r-', 'LineWidth', 2);
    title('Row Norms');
    xlabel('Acquisition Index'); ylabel('Norm');
    grid on;
    
    % Simulated measurements
    subplot(2, 3, 6);
    plot(simulated_measurements, 'g-', 'LineWidth', 2);
    title('Simulated Measurements');
    xlabel('Acquisition Index'); ylabel('Amplitude');
    grid on;
    
    % Save diagnostic plot
    saveas(gcf, fullfile(output_folder, 'h_matrix_diagnostics.png'));
    fprintf('  H matrix diagnostics saved to: h_matrix_diagnostics.png\n');
    
    fprintf('  === END H MATRIX DIAGNOSTICS ===\n');
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

function [reconstructed_image, reconstruction_time, reconstruction_info] = perform_admm_reconstruction_with_diagnostics(measurements, H_matrix, config, params, imaging_grid)
    % Perform ADMM reconstruction with detailed diagnostics
    
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
    
    % ADMM iterations with diagnostics
    primal_residuals = zeros(max_iter, 1);
    dual_residuals = zeros(max_iter, 1);
    objective_values = zeros(max_iter, 1);
    
    fprintf('    Starting ADMM iterations...\n');
    
    for iter = 1:max_iter
        % x-update (least squares with regularization)
        x_prev = x;
        x = (HtH + rho * eye(num_pixels)) \ (Hty + rho * (z - u));
        
        % z-update (TV proximal operator)
        z_prev = z;
        z = prox_tv(x + u, lambda_tv / rho);
        
        % u-update (dual variable)
        u = u + (x - z);
        
        % Calculate residuals and objective
        primal_residual = norm(x - z);
        dual_residual = rho * norm(z - z_prev);
        
        primal_residuals(iter) = primal_residual;
        dual_residuals(iter) = dual_residual;
        
        % Calculate objective value
        data_term = 0.5 * norm(H_matrix * x - measurements)^2;
        tv_term = lambda_tv * tv_norm(x, num_x, num_z);
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
    
    fprintf('    ADMM completed in %.2f seconds\n', reconstruction_time);
    fprintf('    - Iterations: %d\n', iter);
    fprintf('    - Converged: %s\n', mat2str(reconstruction_info.converged));
    fprintf('    - Final primal residual: %.2e\n', primal_residual);
end

function z = prox_tv(x, lambda)
    % TV proximal operator using soft thresholding
    % This is a simplified version - in practice you'd use a proper TV solver
    
    % For now, use soft thresholding on gradients
    [num_x, num_z] = size(reshape(x, sqrt(length(x)), []));
    x_reshaped = reshape(x, num_x, num_z);
    
    % Calculate gradients
    grad_x = diff(x_reshaped, 1, 1);
    grad_z = diff(x_reshaped, 1, 2);
    
    % Soft thresholding
    grad_x_thresh = sign(grad_x) .* max(abs(grad_x) - lambda, 0);
    grad_z_thresh = sign(grad_z) .* max(abs(grad_z) - lambda, 0);
    
    % Reconstruct from gradients (simplified)
    z_reshaped = zeros(size(x_reshaped));
    z_reshaped(1:end-1, :) = z_reshaped(1:end-1, :) + grad_x_thresh;
    z_reshaped(:, 1:end-1) = z_reshaped(:, 1:end-1) + grad_z_thresh;
    
    z = z_reshaped(:);
end

function tv_val = tv_norm(x, num_x, num_z)
    % Calculate TV norm
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

function analyze_admm_results(results, output_folder)
    % Analyze ADMM parameter sweep results
    
    fprintf('  === ADMM RESULTS ANALYSIS ===\n');
    
    % Convert to table for easier analysis
    num_results = length(results);
    lambda_tv = zeros(num_results, 1);
    rho = zeros(num_results, 1);
    max_iter = zeros(num_results, 1);
    psnr = zeros(num_results, 1);
    correlation = zeros(num_results, 1);
    relative_error = zeros(num_results, 1);
    iterations = zeros(num_results, 1);
    converged = zeros(num_results, 1);
    reconstruction_time = zeros(num_results, 1);
    
    for i = 1:num_results
        lambda_tv(i) = results(i).lambda_tv;
        rho(i) = results(i).rho;
        max_iter(i) = results(i).max_iter;
        psnr(i) = results(i).psnr;
        correlation(i) = results(i).correlation;
        relative_error(i) = results(i).relative_error;
        iterations(i) = results(i).iterations;
        converged(i) = results(i).converged;
        reconstruction_time(i) = results(i).reconstruction_time;
    end
    
    % Find best results
    [best_psnr, best_psnr_idx] = max(psnr);
    [best_correlation, best_correlation_idx] = max(correlation);
    [best_converged, best_converged_idx] = max(converged);
    
    fprintf('  Best PSNR: %.2f dB (lambda_tv=%.2f, rho=%.1f, max_iter=%d)\n', ...
        best_psnr, lambda_tv(best_psnr_idx), rho(best_psnr_idx), max_iter(best_psnr_idx));
    fprintf('  Best Correlation: %.4f (lambda_tv=%.2f, rho=%.1f, max_iter=%d)\n', ...
        best_correlation, lambda_tv(best_correlation_idx), rho(best_correlation_idx), max_iter(best_correlation_idx));
    fprintf('  Convergence rate: %.1f%% (%d/%d converged)\n', ...
        mean(converged)*100, sum(converged), num_results);
    
    % Create analysis plots
    figure('Position', [100, 100, 1400, 1000]);
    
    % PSNR vs parameters
    subplot(2, 3, 1);
    scatter3(lambda_tv, rho, psnr, 50, psnr, 'filled');
    colorbar;
    title('PSNR vs Parameters');
    xlabel('Lambda TV'); ylabel('Rho'); zlabel('PSNR (dB)');
    
    % Correlation vs parameters
    subplot(2, 3, 2);
    scatter3(lambda_tv, rho, correlation, 50, correlation, 'filled');
    colorbar;
    title('Correlation vs Parameters');
    xlabel('Lambda TV'); ylabel('Rho'); zlabel('Correlation');
    
    % Convergence analysis
    subplot(2, 3, 3);
    converged_lambda = lambda_tv(converged == 1);
    converged_rho = rho(converged == 1);
    not_converged_lambda = lambda_tv(converged == 0);
    not_converged_rho = rho(converged == 0);
    
    scatter(converged_lambda, converged_rho, 100, 'g', 'filled', 'DisplayName', 'Converged');
    hold on;
    scatter(not_converged_lambda, not_converged_rho, 100, 'r', 'filled', 'DisplayName', 'Not Converged');
    xlabel('Lambda TV'); ylabel('Rho');
    title('Convergence Analysis');
    legend;
    grid on;
    
    % Iterations vs parameters
    subplot(2, 3, 4);
    scatter3(lambda_tv, rho, iterations, 50, iterations, 'filled');
    colorbar;
    title('Iterations vs Parameters');
    xlabel('Lambda TV'); ylabel('Rho'); zlabel('Iterations');
    
    % Reconstruction time vs parameters
    subplot(2, 3, 5);
    scatter3(lambda_tv, rho, reconstruction_time, 50, reconstruction_time, 'filled');
    colorbar;
    title('Reconstruction Time vs Parameters');
    xlabel('Lambda TV'); ylabel('Rho'); zlabel('Time (s)');
    
    % Relative error vs parameters
    subplot(2, 3, 6);
    scatter3(lambda_tv, rho, relative_error, 50, relative_error, 'filled');
    colorbar;
    title('Relative Error vs Parameters');
    xlabel('Lambda TV'); ylabel('Rho'); zlabel('Relative Error');
    
    % Save analysis plot
    saveas(gcf, fullfile(output_folder, 'admm_parameter_analysis.png'));
    
    % Save results to CSV
    results_table = table(lambda_tv, rho, max_iter, psnr, correlation, relative_error, ...
        iterations, converged, reconstruction_time);
    writetable(results_table, fullfile(output_folder, 'admm_results.csv'));
    
    fprintf('  Analysis saved to: admm_parameter_analysis.png\n');
    fprintf('  Results saved to: admm_results.csv\n');
    fprintf('  === END ADMM RESULTS ANALYSIS ===\n');
end

fprintf('\n=== ADMM DIAGNOSTICS COMPLETE ===\n');
fprintf('All results saved to: %s\n', output_folder); 