% =========================================================================
% UNIFIED PMUT ULTRASOUND SIMULATION SCRIPT (v1.10 - Final Optimized)
%
% Description:
% This is the definitive, optimized simulation script. It incorporates all
% successful strategies to produce the highest-quality reconstructions.
%
% v1.10 Improvements:
% - Implements random sparse transmit AND receive apertures for maximum
%   geometric diversity.
% - Re-instates the full test suite of target patterns from v1.5.
% - Uses the memory-safe and CPU-efficient "forward operator" method for
%   reconstruction, eliminating the H-matrix bottleneck.
% =========================================================================
clear; clc; close all;

%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== pMUT SIMULATION (v1.10 - Final Optimized) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('final_simulation_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- Core Physical and Simulation Parameters ---
params = struct();
params.c = 343;
params.fs = 1e6;
params.pmut_width_m = 0.020;
params.kerf_m = 0.005;
params.excitation_amplitude = 1e12;
params.target_SNR_db = 35;

% --- H-Matrix Generation Parameters ---
params.num_acquisitions = 150;
params.num_active_tx = 5; % Number of transmitters to activate
params.num_active_rx = 1; % Number of receivers to activate
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

% --- Target Scene Sweep Parameters ---
target_patterns_to_test = {'3x3_grid', '2x2_grid', 'line_5', 'cross_5'};
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

%% ===== 2. GENERATE H-MATRIX COMPONENTS =====
fprintf('\n--- Generating H-Matrix Components for Forward Operator ---\n');
tic;
[h_components, imaging_grid] = generate_h_matrix_components(params);
fprintf('H-component generation complete. Time: %.2f seconds.\n', toc);

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
    b_vector = run_forward_simulation(h_components, v_true_vector, params);
    
    fprintf('\n--- Reconstructing Image via ADMM ---\n');
    reconstructed_image = run_admm_reconstruction(h_components, b_vector, scene_matrix, params);
    
    fprintf('\n--- Analyzing and Saving Results ---\n');
    analyze_and_plot_results(reconstructed_image, scene_matrix, params, run_output_folder);
end
fprintf('\n\n=== SIMULATION SWEEP COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

function [all_h_aligned, imaging_grid] = generate_h_matrix_components(config)
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    num_x_grid = 4; num_y_grid = 4;
    total_elements = num_x_grid * num_y_grid;
    
    FullAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, ones(num_y_grid,num_x_grid), 1, 1, [0 0 0]);
    
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
        % --- *** NEW: Randomly select Tx AND Rx elements *** ---
        active_tx_indices = randperm(total_elements, config.num_active_tx);
        active_rx_indices = randperm(total_elements, config.num_active_rx);
        
        num_active_tx = length(active_tx_indices);
        
        % Create temporary apertures for this specific shot
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix(active_tx_indices) = 1;
        TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        rx_enabled_matrix(active_rx_indices) = 1;
        RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
        xdc_impulse(RxAperture, ones(1,10));

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
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        
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
    
    disp('  Aligning H-matrix components to common time axis...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), error('No valid H-matrix components were generated.'); end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    
    all_h_aligned = cell(config.num_acquisitions, 1);
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            all_h_aligned{r_acq} = interp1(t_current, all_h_data{r_acq}, t_common_axis, 'linear', 0);
        else
            all_h_aligned{r_acq} = zeros(length(t_common_axis), N_pixels);
        end
    end
    
    xdc_free(FullAperture);
    field_end();
end

function b_vector = run_forward_simulation(h_components, v_true_vector, params)
    num_acqs = params.num_acquisitions;
    num_samples_per_acq = size(h_components{1}, 1);
    Hv_signal = zeros(num_samples_per_acq * num_acqs, 1);
    for i = 1:num_acqs
        start_idx = (i-1) * num_samples_per_acq + 1;
        end_idx = i * num_samples_per_acq;
        Hv_signal(start_idx:end_idx) = h_components{i} * v_true_vector;
    end
    
    signal_power = mean(Hv_signal(:).^2);
    target_SNR_linear = 10^(params.target_SNR_db / 10);
    noise_variance = signal_power / target_SNR_linear;
    noise_sigma = sqrt(noise_variance);
    noise = noise_sigma * randn(size(Hv_signal));
    b_vector = Hv_signal + noise;
end

function reconstructed_image = run_admm_reconstruction(h_components, b_vector, scene_matrix, config)
    imageResolution = size(scene_matrix);
    
    % Find normalization factor from the components
    max_h_val = 0;
    for i = 1:length(h_components)
        max_h_val = max(max_h_val, max(abs(h_components{i}(:))));
    end
    h_norm_factor = max_h_val;
    if h_norm_factor < eps, h_norm_factor = 1; end

    [Afun, Atfun, AtAfun] = operator_setup(h_components, imageResolution, h_norm_factor, config);
    
    b_admm_vec = b_vector(:) / h_norm_factor;
    
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
    
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    fprintf('  Starting ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 5) == 0, fprintf('  ADMM iteration %d/%d...\n', k_admm, config.admm_max_iter); end
        x_prev = x_admm_img_iter;
        
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
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

function [Afun, Atfun, AtAfun] = operator_setup(h_components, imageResolution, h_norm_factor, config)
    num_acqs = config.num_acquisitions;
    num_pixels = prod(imageResolution);
    num_samples_per_acq = size(h_components{1}, 1);

    function y = forward_op(x_img)
        x_vec = x_img(:);
        y = zeros(num_samples_per_acq * num_acqs, 1);
        for i = 1:num_acqs
            start_idx = (i-1) * num_samples_per_acq + 1;
            end_idx = i * num_samples_per_acq;
            y(start_idx:end_idx) = (h_components{i} * x_vec) / h_norm_factor;
        end
    end

    function x_img = transpose_op(y_vec)
        y_reshaped = reshape(y_vec, [num_samples_per_acq, num_acqs]);
        x_vec = zeros(num_pixels, 1);
        for i = 1:num_acqs
            x_vec = x_vec + (h_components{i}' * y_reshaped(:, i)) / h_norm_factor;
        end
        x_img = reshape(x_vec, imageResolution);
    end

    Afun = @(x) forward_op(x);
    Atfun = @(y) transpose_op(y);
    AtAfun = @(x) transpose_op(forward_op(x));
end

% (Other helper functions like create_scene, difference_operators, etc. are unchanged)
% ...
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

function analyze_and_plot_results(reconstructed_image, scene_matrix, params, output_folder)
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(abs(reconstructed_image(:)));
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    fprintf('  PSNR: %.2f dB\n', psnr);
    
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
    
    metrics = struct('psnr', psnr);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
end
