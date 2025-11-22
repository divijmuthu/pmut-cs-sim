% =========================================================================
% UNIFIED PMUT ULTRASOUND SIMULATION SCRIPT (v1.27 - Decoupled Antipodal Rx)
%
% Description:
% This version implements a robust, intelligent setup for the antipodal
% receiver strategy. To prevent selection conflicts, it first reserves a
% random antipodal pair for receiving, then selects the 8 transmitters
% from the remaining 14 elements using the Balanced Quadrant Pairs method.
% =========================================================================
clear; clc; close all;

% BEST CASE IS MORE SPACE BETWEEN PMUTS IT REALLY REALLY HELPS
% 2.5 cm gap was excellent

%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== pMUT SIMULATION (v1.27 - Decoupled Antipodal Rx) ===\n');
% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('decoupled_antipodal_rx_sim_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);
% --- ADMM HYPERPARAMETER TUNING ---
params.lambda_tv_reg = 5.0;
params.rho_admm = 6.73;
% --- Core Physical and Simulation Parameters ---
params.c = 343;
params.fs = 1e6;
params.pmut_width_m = 0.002;
params.kerf_m = 0.015; % real thing is 0.005, 0.025 is excellent, 
% 0.01 or 0.015 are prob more realistic esp for current pcb
% 0.01 gives >20 dB for first case which is meh 
% 0.015 has a clear improvement to coherence prob reconst too
    % 27 dB SNR which is quite nice maybe 0.0125 is good
    % discuss with NL more ofc also unoptimized ADMM too
    % probably too big for real pcb tho 2mm 
params.excitation_amplitude = 1e12;
params.target_SNR_db = 35;
% --- H-Matrix Generation Parameters ---
params.num_acquisitions = 150;
params.num_active_tx = 8;
params.num_active_rx = 2;
params.max_delay_rand_us = 500;
params.apodization_mode = 'random';
params.f_min_hz = 45e3;
params.f_max_hz = 65e3;
params.use_phase_shifts = true;
params.assembly_chunk_size = 25;
% --- Imaging Grid Parameters ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;
% --- Target Scene Sweep Parameters ---
target_patterns_to_test = {'3x3_grid', '2x2_grid', 'line_5', 'cross_5'};
params.grid_spacing_mm = 25;
params.target_sizes_mm = [3, 3, 4, 4, 4, 5, 5, 5, 6];
params.target_intensities = [0.8, 0.9, 1.0, 1.0, 1.1, 1.1, 1.2, 1.2, 1.3];
% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;
%% ===== 2. GENERATE H-MATRIX AND ANALYZE COHERENCE =====
fprintf('\n--- Generating Full H-Matrix ---\n');
tic;
[H, imaging_grid] = generate_h_matrix(params);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);
fprintf('\n--- Analyzing H-Matrix Coherence ---\n');
[max_coherence, coherence_matrix] = analyze_coherence(H);
fprintf('  Max Coherence: %.4f\n', max_coherence);
figure('Visible', 'on', 'Position', [100, 100, 700, 600]);
imagesc(coherence_matrix); axis square; colorbar;
title(sprintf('H-Matrix Coherence\nMax: %.4f', max_coherence));
xlabel('Column Index'); ylabel('Column Index');
saveas(gcf, fullfile(output_folder, 'coherence_analysis.png'));
%% ===== 3. SEQUENTIAL SCENE SIMULATION SWEEP =====
for i = 1:length(target_patterns_to_test)
    current_pattern = target_patterns_to_test{i};
    fprintf('\n============================================================\n');
    fprintf('--- Running Simulation for Scene %d/%d: "%s" ---\n', i, length(target_patterns_to_test), current_pattern);
    fprintf('============================================================\n');
    
    run_name = sprintf('run_%02d_%s', i, current_pattern);
    run_output_folder = fullfile(output_folder, run_name);
    if ~exist(run_output_folder, 'dir'), mkdir(run_output_folder); end
    
    params.target_pattern = current_pattern;
    
    fprintf('\n--- Creating Target Scene ---\n');
    scene_matrix = create_scene(imaging_grid, params);
    v_true_vector = scene_matrix(:);
    visualize_scene(scene_matrix, imaging_grid, run_output_folder);
    
    fprintf('\n--- Running Forward Simulation ---\n');
    b_vector = run_forward_simulation(H, v_true_vector, params);
    
    fprintf('\n--- Reconstructing Image via ADMM ---\n');
    reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, params);
    
    fprintf('\n--- Analyzing and Saving Results ---\n');
    analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, run_output_folder);
end
fprintf('\n\n=== SIMULATION SWEEP COMPLETE ===\n');
%% ===== HELPER FUNCTIONS =====
function [H, imaging_grid] = generate_h_matrix(config)
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    num_x_grid = 4; num_y_grid = 4;
    total_elements = num_x_grid * num_y_grid;
    
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
        
        % ================================================================
        % --- NEW: Decoupled Antipodal Selection Logic (v1.27) ---
        
        % 1. Define all possible antipodal pairs. The opposite of index i is 17-i.
        all_rx_pairs = [ (1:8)', (17-(1:8))' ];
        
        % 2. FIRST, select a random pair to be the receivers. This reserves them.
        chosen_pair_idx = randi(size(all_rx_pairs, 1));
        active_rx_indices = all_rx_pairs(chosen_pair_idx, :);
        
        % 3. THEN, select transmitters from the remaining 14 elements.
        available_for_tx = setdiff(1:total_elements, active_rx_indices);
        
        % 4. Dynamically define quadrants based on available TX elements
        quad_defs = {[1, 2, 5, 6], [3, 4, 7, 8], [9, 10, 13, 14], [11, 12, 15, 16]};
        quadrants_available = cell(1,4);
        for q = 1:4
            quadrants_available{q} = intersect(quad_defs{q}, available_for_tx);
        end
        
        % 5. Run the Balanced Quadrant Pairs logic on the available elements
        active_tx_indices = zeros(1, config.num_active_tx);
        elements_per_quad = config.num_active_tx / 4;
        currentIndex = 1;
        for q = 1:4
            quad_elements = quadrants_available{q};
            if length(quad_elements) < elements_per_quad
                 error('A quadrant does not have enough available elements for transmitter selection after reserving receivers.');
            end
            rand_indices = randperm(length(quad_elements), elements_per_quad);
            selected_elements = quad_elements(rand_indices);
            active_tx_indices(currentIndex : currentIndex + elements_per_quad - 1) = selected_elements;
            currentIndex = currentIndex + elements_per_quad;
        end
        % ================================================================
        
        num_active_tx = length(active_tx_indices);
        
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix(active_tx_indices) = 1;
        TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        rx_enabled_matrix(active_rx_indices) = 1;
        RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
        
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        
        f0 = mean(tx_frequencies);
        t_pulse = 0:1/fs:(5/f0);
        impulse_response = sin(2*pi*f0*t_pulse) .* tukeywin(length(t_pulse), 0.4)';
        xdc_impulse(TxAperture, impulse_response * config.excitation_amplitude);
        xdc_impulse(RxAperture, ones(1,10));
        
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        
        delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        
        xdc_apodization(TxAperture, 0, apod_weights);
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        xdc_free(TxAperture); xdc_free(RxAperture);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, config, N_pixels);
    
    field_end();
end
function H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, config, N_pixels)
    disp('  Assembling H-matrix using chunked triplet method...');
    fs = config.fs;
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), H = sparse(config.num_acquisitions, N_pixels); return; end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    
    chunk_size = config.assembly_chunk_size;
    num_chunks = ceil(config.num_acquisitions / chunk_size);
    H_chunks = cell(num_chunks, 1);
    
    wb_chunk = waitbar(0, 'Assembling H-matrix in chunks...');
    for c_idx = 1:num_chunks
        start_acq = (c_idx - 1) * chunk_size + 1;
        end_acq = min(c_idx * chunk_size, config.num_acquisitions);
        num_acqs_in_chunk = end_acq - start_acq + 1;
        
        total_nnz_chunk = 0;
        all_h_aligned_chunk = cell(num_acqs_in_chunk, 1);
        for r_acq_local = 1:num_acqs_in_chunk
            r_acq_global = start_acq + r_acq_local - 1;
            if all_K_values(r_acq_global) > 0 && ~isempty(all_h_data{r_acq_global})
                t_current = all_start_times(r_acq_global) + (0:(all_K_values(r_acq_global) - 1)) / fs;
                h_aligned = interp1(t_current, all_h_data{r_acq_global}, t_common_axis, 'linear', 0);
                total_nnz_chunk = total_nnz_chunk + nnz(h_aligned);
                all_h_aligned_chunk{r_acq_local} = h_aligned;
            end
        end
        
        I_list = zeros(total_nnz_chunk, 1); J_list = zeros(total_nnz_chunk, 1); S_list = zeros(total_nnz_chunk, 1);
        currentIndex = 1;
        for r_acq_local = 1:num_acqs_in_chunk
            if ~isempty(all_h_aligned_chunk{r_acq_local})
                [row_idx, col_idx, s_vals] = find(all_h_aligned_chunk{r_acq_local});
                global_row_idx = row_idx + (r_acq_local - 1) * K_global_per_acq;
                num_elements = length(s_vals);
                I_list(currentIndex : currentIndex + num_elements - 1) = global_row_idx;
                J_list(currentIndex : currentIndex + num_elements - 1) = col_idx;
                S_list(currentIndex : currentIndex + num_elements - 1) = s_vals;
                currentIndex = currentIndex + num_elements;
            end
        end
        
        total_rows_chunk = K_global_per_acq * num_acqs_in_chunk;
        H_chunks{c_idx} = sparse(I_list, J_list, S_list, total_rows_chunk, N_pixels);
        waitbar(c_idx / num_chunks, wb_chunk, sprintf('Assembled chunk %d/%d', c_idx, num_chunks));
    end
    close(wb_chunk);
    
    H = vertcat(H_chunks{:});
end
function [max_coherence, coherence_matrix] = analyze_coherence(H)
    fprintf('  Calculating coherence (this may be slow)...');
    tic;
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0; coherence_matrix = zeros(size(H,2));
    end
    fprintf(' Done. (%.2f seconds)\n', toc);
end
function scene_matrix = create_scene(imaging_grid, params)
    scene_matrix = zeros(size(imaging_grid.X_mesh));
    X_mm = imaging_grid.X_mesh * 1000;
    Y_mm = imaging_grid.Y_mesh * 1000;
    
    grid_spacing_mm = params.grid_spacing_mm;
    
    switch params.target_pattern
        case '3x3_grid'
            positions = [];
            for row = -1:1, for col = -1:1, positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm]; end, end
        case '2x2_grid'
            positions = [];
            for row = -0.5:1:0.5, for col = -0.5:1:0.5, positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm]; end, end
        case 'line_5'
            positions = [(-2:2)' * grid_spacing_mm, zeros(5,1)];
        case 'cross_5'
            positions = [0,0; -1,0; 1,0; 0,-1; 0,1] * grid_spacing_mm;
    end
    
    num_targets = size(positions, 1);
    target_sizes = repmat(params.target_sizes_mm, 1, ceil(num_targets/length(params.target_sizes_mm)));
    target_intensities = repmat(params.target_intensities, 1, ceil(num_targets/length(params.target_intensities)));
    
    for i = 1:num_targets
        x_pos_mm = positions(i, 1); y_pos_mm = positions(i, 2);
        target_size_mm = target_sizes(i); target_intensity = target_intensities(i);
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iy_scene] = min(abs(Y_mm(:,1) - y_pos_mm));
        target_radius_pixels = round(target_size_mm / (2 * (params.grid_step_m*1000)));
        if target_radius_pixels < 1, target_radius_pixels = 1; end
        for dy = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                ix_target = ix_scene + dx; iy_target = iy_scene + dy;
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

        if ishandle(wb_admm)
            waitbar(k_admm / config.admm_max_iter, wb_admm);
        end
    end

    if ishandle(wb_admm)
        close(wb_admm);
    end

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
    
    [max_coherence, ~] = analyze_coherence(H);
    
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
    sgtitle(sprintf('Decoupled Antipodal Rx (Tx=%d, Rx=%d, lambda=%.1f)', params.num_active_tx, params.num_active_rx, params.lambda_tv_reg));
    saveas(gcf, fullfile(output_folder, 'reconstruction_comparison.png'));
    
    metrics = struct('psnr', psnr, 'max_coherence', max_coherence, 'params', params);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
end