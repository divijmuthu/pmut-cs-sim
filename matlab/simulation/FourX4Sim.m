% =========================================================================
% UNIFIED PMUT ULTRASOUND SIMULATION SCRIPT (v1.9 - 4x4 Grid)
%
% Description:
% This version adapts the successful v1.5 simulation framework to test a
% physically realistic hardware setup. It simulates a fixed 4x4 grid of
% 16 pMUTs and, for each acquisition, randomly activates a sparse subset
% of 5 elements. This directly tests if a practical array can provide the
% geometric diversity needed for a low-coherence H-matrix and a high-
% quality reconstruction.
%
% v1.9 Improvements:
% - Replaced the 10x10 "virtual grid" with a fixed 4x4 physical grid.
% - Randomly selects 5 out of the 16 available elements for each shot.
% =========================================================================
clear; clc; close all;

%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== pMUT SIMULATION (v1.9 - 4x4 Grid) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('4x4_grid_simulation_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- Core Physical and Simulation Parameters ---
params = struct();
params.c = 343;
params.fs = 1e6;
params.pmut_width_m = 0.020;
params.kerf_m = 0.005; % Kerf between elements in the 4x4 grid
params.excitation_amplitude = 1e12;
params.target_SNR_db = 35;

% --- H-Matrix Generation Parameters ---
params.num_acquisitions = 150;
params.num_active_tx = 5; % Number of elements to activate from the 4x4 grid
params.max_delay_rand_us = 500;
params.apodization_mode = 'random';
params.f_min_hz = 45e3;
params.f_max_hz = 65e3;
params.use_phase_shifts = true;

% --- Imaging Grid Parameters ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;

% --- Target Scene Parameters ---
params.target_pattern = '3x3_grid';
params.grid_spacing_mm = 25;
params.target_sizes_mm = [3, 3, 4, 4, 4, 5, 5, 5, 6];
params.target_intensities = [0.8, 0.9, 1.0, 1.0, 1.1, 1.1, 1.2, 1.2, 1.3];

% --- ADMM Reconstruction Parameters (from v1.5) ---
params.rho_admm = 6.73;
params.lambda_tv_reg = 1.2;
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;

%% ===== 2. GENERATE H-MATRIX =====
fprintf('\n--- Generating H-Matrix with 4x4 Random Subset Selection ---\n');
tic;
[H, imaging_grid] = generate_h_matrix(params);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

%% ===== 3. CREATE TARGET SCENE =====
fprintf('\n--- Creating Target Scene ---\n');
tic;
scene_matrix = create_scene(imaging_grid, params);
v_true_vector = scene_matrix(:);
fprintf('Scene creation complete. Time: %.2f seconds.\n', toc);
visualize_scene(scene_matrix, imaging_grid, output_folder);

%% ===== 4. RUN FORWARD SIMULATION (b = H*v + n) =====
fprintf('\n--- Running Forward Simulation ---\n');
tic;
b_vector = run_forward_simulation(H, v_true_vector, params);
fprintf('Forward simulation complete. Time: %.2f seconds.\n', toc);

%% ===== 5. RECONSTRUCT IMAGE USING ADMM =====
fprintf('\n--- Reconstructing Image via ADMM ---\n');
tic;
reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);

%% ===== 6. ANALYZE AND SAVE RESULTS =====
fprintf('\n--- Analyzing and Saving Results ---\n');
analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, output_folder);
fprintf('\n\n=== SIMULATION COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix(config)
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    % --- *** NEW: Define a fixed 4x4 physical array *** ---
    num_x_grid = 4;
    num_y_grid = 4;
    total_elements = num_x_grid * num_y_grid;
    
    % Define the full aperture containing all 16 possible pMUT locations
    FullTxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, ones(num_y_grid,num_x_grid), 1, 1, [0 0 0]);
    
    % Define a single receiver at the center of the coordinate system
    RxAperture = xdc_2d_array(1, 1, config.pmut_width_m, config.pmut_width_m, 0, 0, 1, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
    y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        % --- *** NEW: Randomly select 5 out of the 16 available elements *** ---
        active_indices = randperm(total_elements, config.num_active_tx);
        num_active_tx = length(active_indices);
        
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        
        tx_signals = cell(1, num_active_tx); max_len = 0;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k); t = 0:1/fs:duration;
            random_phase = 0;
            if config.use_phase_shifts, random_phase = 2 * pi * rand(); end
            signal_base = sin(2 * pi * tx_frequencies(k) * t + random_phase);
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx, sig = tx_signals{k}; composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig; end
        
        xdc_impulse(FullTxAperture, composite_waveform * config.excitation_amplitude);
        
        % Create vectors that are the full size of the 4x4 grid
        full_apod_vector = zeros(1, total_elements);
        full_delay_vector = zeros(1, total_elements);
        
        % Generate random apodization weights for the active elements
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        
        % Place the parameters for the active elements into the full vectors
        full_apod_vector(active_indices) = apod_weights;
        full_delay_vector(active_indices) = delays_us * 1e-6;
        
        xdc_apodization(FullTxAperture, 0, full_apod_vector);
        xdc_focus_times(FullTxAperture, 0, full_delay_vector);
        
        [h_r, start_time_r] = calc_hhp(FullTxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    H = assemble_h_matrix(all_h_data, all_start_times, all_K_values, config, N_pixels);
    
    xdc_free(FullTxAperture); xdc_free(RxAperture); field_end();
end

function H = assemble_h_matrix(all_h_data, all_start_times, all_K_values, config, N_pixels)
    disp('  Assembling H-matrix using interpolation...');
    fs = config.fs;
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), H = sparse(config.num_acquisitions, N_pixels); return; end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0, K_global_per_acq = 1; end
    total_rows = K_global_per_acq * config.num_acquisitions;
    H = spalloc(total_rows, N_pixels, round(sum(all_K_values) * N_pixels * 0.05));
    current_row_offset = 0;
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned = interp1(t_current, all_h_data{r_acq}, t_common_axis, 'linear', 0);
            row_indices = current_row_offset + (1:K_global_per_acq);
            if max(row_indices) <= size(H, 1)
                H(row_indices, :) = h_aligned;
            end
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
end

function scene_matrix = create_scene(imaging_grid, params)
    scene_matrix = zeros(size(imaging_grid.X_mesh));
    X_mm = imaging_grid.X_mesh * 1000;
    Y_mm = imaging_grid.Y_mesh * 1000;
    
    grid_spacing_mm = params.grid_spacing_mm;
    
    switch params.target_pattern
        case '3x3_grid'
            positions = [];
            for row = -1:1, for col = -1:1
                positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm];
            end, end
    end
    
    num_targets = size(positions, 1);
    target_sizes = repmat(params.target_sizes_mm, 1, ceil(num_targets/length(params.target_sizes_mm)));
    target_intensities = repmat(params.target_intensities, 1, ceil(num_targets/length(params.target_intensities)));
    
    for i = 1:num_targets
        x_pos_mm = positions(i, 1);
        y_pos_mm = positions(i, 2);
        
        target_size_mm = target_sizes(i);
        target_intensity = target_intensities(i);
        
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iy_scene] = min(abs(Y_mm(:,1) - y_pos_mm));
        
        target_radius_pixels = round(target_size_mm / (2 * (params.grid_step_m*1000)));
        if target_radius_pixels < 1, target_radius_pixels = 1; end
        
        for dy = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                ix_target = ix_scene + dx;
                iy_target = iy_scene + dy;
                if ix_target > 0 && ix_target <= size(X_mm, 2) && iy_target > 0 && iy_target <= size(Y_mm, 1)
                    scene_matrix(iy_target, ix_target) = target_intensity;
                end
            end
        end
    end
end

function b_vector = run_forward_simulation(H, v_true_vector, params)
    Hv_signal = H * v_true_vector;
    signal_power = mean(Hv_signal(:).^2);
    target_SNR_linear = 10^(params.target_SNR_db / 10);
    noise_variance = signal_power / target_SNR_linear;
    noise_sigma = sqrt(noise_variance);
    noise = noise_sigma * randn(size(Hv_signal));
    b_vector = Hv_signal + noise;
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config)
    imageResolution = size(scene_matrix);
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor;
    At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    fprintf('  Starting ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 5) == 0, fprintf('  ADMM iteration %d/%d...\n', k_admm, config.admm_max_iter); end
        x_prev = x_admm_img_iter;
        
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun_admm_img(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
        [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter);
        x_admm_img_iter = reshape(x_vec_new, imageResolution);
        
        kap = config.lambda_tv_reg / config.rho_admm;
        v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
        v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm = max(v_norm, eps);
        shr = max(0, 1 - kap ./ v_norm);
        z_admm_grad_iter = v_z_upd .* shr;
        u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
        
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol, fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        end
        waitbar(k_admm / config.admm_max_iter, wb_admm);
    end
    close(wb_admm);
    reconstructed_image = x_admm_img_iter;
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    Afun_admm = @(x) A_admm * x(:); Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [Dx, Dy] = difference_operators(imageSize)
   rows = imageSize(1); cols = imageSize(2); N_img_pixels = rows * cols;
   Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
   Dx( (cols-1)*rows+1 : cols*rows , :) = 0;
   Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
   Dy( rows:rows:N_img_pixels , :) = 0;
end

function visualize_scene(scene_matrix, imaging_grid, output_folder)
    figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
    imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Y_mesh(:,1)*1000, scene_matrix);
    colormap(gray); colorbar; axis image;
    title('Ground Truth Target Scene');
    xlabel('X Position (mm)'); ylabel('Y Position (mm)');
    saveas(gcf, fullfile(output_folder, 'ground_truth_scene.png'));
    close(gcf);
end

function analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, output_folder)
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(abs(reconstructed_image(:)));
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0;
    end
    
    fprintf('  PSNR: %.2f dB\n', psnr);
    fprintf('  Max Coherence: %.4f\n', max_coherence);
    
    figure('Visible', 'on', 'Position', [100, 100, 1200, 500]);
    subplot(1, 2, 1);
    imagesc(scene_matrix); colormap(gray); axis image; colorbar;
    title('Ground Truth');
    xlabel('X Index'); ylabel('Y Index');
    subplot(1, 2, 2);
    imagesc(reconstructed_image); colormap(gray); axis image; colorbar;
    title(sprintf('Reconstruction (PSNR: %.2f dB)', psnr));
    xlabel('X Index'); ylabel('Y Index');
    sgtitle('Simulation Results');
    saveas(gcf, fullfile(output_folder, 'reconstruction_comparison.png'));
    
    metrics = struct('psnr', psnr, 'max_coherence', max_coherence);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
end
