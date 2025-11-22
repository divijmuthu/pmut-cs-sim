%% test_fixed_parameter_sweep.m - Test the fixed parameter sweep functions
% Test the core functions without running the full sweep

clear; clc; close all;

fprintf('=== TESTING FIXED PARAMETER SWEEP ===\n');

% Test basic configuration
config = struct();
config.c = 343;
config.fs = 1e6;
config.pmut_width_m = 0.020;
config.tx_pool_width_m = 0.200;
config.grid_width_m = 0.150;
config.target_distance_m = 0.150;
config.grid_depth_range_m = 0.100;
config.grid_step_m = 0.0075;  % 7.5mm
config.excitation_amplitude = 1e15;
config.pmut_resonance_freq = 57700;
config.pmut_bandwidth = 2520;
config.impulse_duration_us = 10;
config.num_active_tx = 5;
config.max_delay_rand_us = 500;
config.apodization_mode = 'uniform';
config.frequency_offset_hz = 0;
config.num_acquisitions = 20;
config.numItersADMM = 50;
config.rho_admm = 6.73;
config.lambda_tv_reg = 0.7438;
config.admm_tol = 1.2e-5;
config.admm_max_iter = 50;
config.pcg_max_iter = 30;
config.pcg_tol = 1e-8;
config.target_SNR_db = 35;

% Test parameters
params = struct();
params.grid_step_mm = 7.5;
params.target_size_mm = 4;
params.grid_spacing_mm = 20;
params.num_acquisitions = 20;
params.lambda_tv_reg = 0.7438;
params.target_pattern = '3x3_grid';

fprintf('Testing scene creation...\n');
try
    % Create imaging grid
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    
    scene_matrix = create_scene_with_pattern(X_mesh, Z_mesh, config, params);
    fprintf('✓ Scene creation successful\n');
    fprintf('  Scene size: %d x %d\n', size(scene_matrix));
    fprintf('  Non-zero pixels: %d\n', nnz(scene_matrix));
catch ME
    fprintf('✗ Scene creation failed: %s\n', ME.message);
    return;
end

fprintf('\nTesting synthetic H matrix creation...\n');
try
    % Create a synthetic H matrix for testing (without Field II)
    N_pixels = numel(scene_matrix);
    N_measurements = config.num_acquisitions * 100; % Synthetic measurements
    
    % Create a random sparse H matrix
    H = sprandn(N_measurements, N_pixels, 0.1);
    
    % Normalize columns
    col_norms = sqrt(sum(H.^2, 1));
    H = H ./ col_norms;
    
    fprintf('✓ Synthetic H matrix creation successful\n');
    fprintf('  H matrix size: %d x %d\n', size(H));
    fprintf('  Non-zero elements: %d\n', nnz(H));
catch ME
    fprintf('✗ Synthetic H matrix creation failed: %s\n', ME.message);
    return;
end

fprintf('\nTesting forward simulation...\n');
try
    v_true_vector = scene_matrix(:);
    Hv_signal = H * v_true_vector;
    signal_power = mean(Hv_signal(:).^2);
    target_SNR_linear = 10^(config.target_SNR_db / 10);
    noise_variance = signal_power / target_SNR_linear;
    noise_sigma = sqrt(noise_variance);
    noise = noise_sigma * randn(size(Hv_signal));
    b_vector = Hv_signal + noise;
    fprintf('✓ Forward simulation successful\n');
    fprintf('  Signal power: %.2e\n', signal_power);
    fprintf('  Noise variance: %.2e\n', noise_variance);
catch ME
    fprintf('✗ Forward simulation failed: %s\n', ME.message);
    return;
end

fprintf('\nTesting ADMM reconstruction...\n');
try
    output_folder = 'test_fixed_output';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, output_folder);
    fprintf('✓ ADMM reconstruction successful\n');
    fprintf('  Reconstruction size: %d x %d\n', size(reconstructed_image));
catch ME
    fprintf('✗ ADMM reconstruction failed: %s\n', ME.message);
    fprintf('  Error details: %s\n', getReport(ME, 'extended'));
    return;
end

fprintf('\nTesting metrics calculation...\n');
try
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    scene_vec = scene_norm(:);
    recon_vec = recon_norm(:);
    correlation = (scene_vec - mean(scene_vec))' * (recon_vec - mean(recon_vec)) / ...
        (sqrt(sum((scene_vec - mean(scene_vec)).^2)) * sqrt(sum((recon_vec - mean(recon_vec)).^2)));
    
    fprintf('✓ Metrics calculation successful\n');
    fprintf('  PSNR: %.2f dB\n', psnr);
    fprintf('  Correlation: %.4f\n', correlation);
catch ME
    fprintf('✗ Metrics calculation failed: %s\n', ME.message);
    return;
end

fprintf('\n=== FIXED TEST COMPLETE ===\n');
fprintf('All components working correctly!\n');

%% Helper Functions (copied from fixed parameter sweep)

function scene_matrix = create_scene_with_pattern(X_mesh, Z_mesh, config, params)
    % Create scene based on specified pattern
    
    scene_matrix = zeros(size(X_mesh));
    X_mm = X_mesh * 1000;
    Z_mm = Z_mesh * 1000;
    
    grid_spacing_mm = params.grid_spacing_mm;
    target_size_mm = params.target_size_mm;
    grid_offset_x_mm = 0;
    grid_offset_z_mm = 150;
    
    grid_step_mm = config.grid_step_m * 1000;
    target_radius_pixels = round(target_size_mm / (2 * grid_step_mm));
    if target_radius_pixels < 1
        target_radius_pixels = 1;
    end
    
    switch params.target_pattern
        case '3x3_grid'
            % 3x3 grid pattern
            target_positions = [];
            for row = 1:3
                for col = 1:3
                    x_pos_mm = grid_offset_x_mm + (col - 2) * grid_spacing_mm;
                    z_pos_mm = grid_offset_z_mm + (row - 2) * grid_spacing_mm;
                    target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
                end
            end
            
        case '2x2_grid'
            % 2x2 grid pattern
            target_positions = [];
            for row = 1:2
                for col = 1:2
                    x_pos_mm = grid_offset_x_mm + (col - 1.5) * grid_spacing_mm;
                    z_pos_mm = grid_offset_z_mm + (row - 1.5) * grid_spacing_mm;
                    target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
                end
            end
            
        case 'line_5'
            % Line of 5 targets
            target_positions = [];
            for i = 1:5
                x_pos_mm = grid_offset_x_mm + (i - 3) * grid_spacing_mm;
                z_pos_mm = grid_offset_z_mm;
                target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
            end
            
        case 'cross_5'
            % Cross pattern with 5 targets
            target_positions = [
                grid_offset_x_mm, grid_offset_z_mm, 1.0;  % Center
                grid_offset_x_mm - grid_spacing_mm, grid_offset_z_mm, 1.0;  % Left
                grid_offset_x_mm + grid_spacing_mm, grid_offset_z_mm, 1.0;  % Right
                grid_offset_x_mm, grid_offset_z_mm - grid_spacing_mm, 1.0;  % Top
                grid_offset_x_mm, grid_offset_z_mm + grid_spacing_mm, 1.0;  % Bottom
            ];
    end
    
    % Create targets
    for i = 1:size(target_positions, 1)
        x_pos_mm = target_positions(i, 1);
        z_pos_mm = target_positions(i, 2);
        amplitude = target_positions(i, 3);
        
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iz_scene] = min(abs(Z_mm(:,1) - z_pos_mm));
        
        for dz = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                ix_target = ix_scene + dx;
                iz_target = iz_scene + dz;
                
                if ix_target > 0 && ix_target <= size(X_mesh, 2) && ...
                   iz_target > 0 && iz_target <= size(X_mesh, 1)
                    scene_matrix(iz_target, ix_target) = amplitude;
                end
            end
        end
    end
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, output_folder)
    % Run ADMM reconstruction with TV regularization (EXACT COPY FROM WORKING SCRIPT)
    
    fprintf('  Starting ADMM reconstruction...\n');
    
    % Setup ADMM problem
    imageResolution = size(scene_matrix);
    N_pixels = numel(scene_matrix);
    
    % Normalize true scene for metrics
    v_true_vec_norm = scene_matrix(:) / max(scene_matrix(:));
    
    % Setup operators (EXACT COPY)
    A_matrix = H;
    if size(A_matrix, 1) ~= length(b_vector)
        fprintf('  ERROR: H matrix rows (%d) != measurement vector length (%d)\n', size(A_matrix, 1), length(b_vector));
        reconstructed_image = zeros(imageResolution);
        return;
    end
    
    % Normalization
    H_norm_factor = max(abs(A_matrix(:)));
    if H_norm_factor < eps
        H_norm_factor = 1;
    end
    A_admm = A_matrix ./ H_norm_factor;
    At_admm = transpose(A_admm);
    b_admm_vec = b_vector(:) / H_norm_factor;
    
    % Operator setup
    [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    
    % Pre-allocate variables (EXACT COPY)
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    % PCG function (EXACT COPY)
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    % Pre-allocate tracking arrays
    PSNR_admm_iters = zeros([config.admm_max_iter 1]);
    residuals_admm_iters = zeros([config.admm_max_iter 1]);
    
    % ADMM iterations (EXACT COPY)
    converged = false;
    
    for k_admm = 1:config.admm_max_iter
        % Store previous iteration
        x_prev = x_admm_img_iter;
        
        % ADMM update
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = admm_iteration(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            config.rho_admm, config.lambda_tv_reg, Hfun_pcg_admm, config);
        
        % Metrics calculation
        [PSNR_admm_iters(k_admm), residuals_admm_iters(k_admm)] = calculate_metrics(...
            x_admm_img_iter, v_true_vec_norm, b_admm_vec, Afun_admm, opDx_tv, config.lambda_tv_reg);
        
        % Check convergence
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol
                converged = true;
                fprintf('    ADMM converged at iteration %d (rel change: %.2e)\n', k_admm, rel_change);
                break;
            end
        end
        
        % Print progress
        if mod(k_admm, 10) == 0
            fprintf('    Iteration %d/%d: PSNR=%.2f dB, Residual=%.2e\n', ...
                k_admm, config.admm_max_iter, PSNR_admm_iters(k_admm), residuals_admm_iters(k_admm));
        end
    end
    
    if ~converged
        fprintf('    ADMM did not converge within %d iterations\n', config.admm_max_iter);
    end
    
    % Final result
    reconstructed_image = x_admm_img_iter;
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    % Setup operators for ADMM
    
    % Matrix operators
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(b) reshape(At_admm * b, imageResolution);
    AtAfun_admm_img = @(x) reshape(At_admm * (A_admm * x(:)), imageResolution);
    
    % TV operators - FIXED DIMENSIONS
    [Dx, Dy] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx * x(:), Dy * x(:)];
    opDtx_tv = @(v) reshape(Dx' * v(:,1) + Dy' * v(:,2), imageResolution);
    opDtDx_tv = @(x) reshape(Dx' * (Dx * x(:)) + Dy' * (Dy * x(:)), imageResolution);
end

function [Dx, Dy] = difference_operators(imageSize)
    % Create difference operators for TV regularization
    
    N = prod(imageSize);
    
    % Forward differences
    Dx = sparse(N, N);
    Dy = sparse(N, N);
    
    for i = 1:imageSize(1)
        for j = 1:imageSize(2)
            idx = sub2ind(imageSize, i, j);
            
            % X difference
            if j < imageSize(2)
                idx_next = sub2ind(imageSize, i, j+1);
                Dx(idx, idx) = -1;
                Dx(idx, idx_next) = 1;
            end
            
            % Y difference
            if i < imageSize(1)
                idx_next = sub2ind(imageSize, i+1, j);
                Dy(idx, idx) = -1;
                Dy(idx, idx_next) = 1;
            end
        end
    end
end

function [x_new, z_new, u_new] = admm_iteration(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, Hfun_pcg_admm, config)
    % Single ADMM iteration (EXACT COPY FROM WORKING SCRIPT)
    
    % x-update
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    
    % z-update
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm = max(v_norm, eps);
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    
    % u-update
    u_new = u_old + opDx_tv(x_new) - z_new;
end

function [PSNR, residual] = calculate_metrics(x_admm_img_iter, v_true_vec_norm, b_admm_vec, Afun_admm, opDx_tv, lambda_tv_reg)
    % Calculate reconstruction metrics
    
    % Normalization
    x_scl = real(x_admm_img_iter(:));
    x_range = max(x_scl) - min(x_scl);
    if x_range > eps
        x_scl = (x_scl - min(x_scl)) / x_range;
    else
        x_scl = zeros(size(x_scl));
    end
    
    % PSNR calculation
    MSE_curr = mean((x_scl - v_true_vec_norm).^2);
    PSNR = 10 * log10(1 / MSE_curr);
    
    % Residual calculation
    r1 = b_admm_vec - Afun_admm(x_admm_img_iter);
    r2 = opDx_tv(x_admm_img_iter);
    tv_n = sum(sqrt(sum(r2.^2, 2)));
    residual = 0.5 * sum(r1(:).^2) + lambda_tv_reg * tv_n;
end 