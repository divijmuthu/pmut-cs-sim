% =========================================================================
% UNIFIED PMUT ULTRASOUND SIMULATION SCRIPT (v1.1)
%
% Description:
% This script integrates the best components from previous development into
% a single, cohesive, end-to-end simulation framework.
%
% v1.1 Improvements:
% - Increased geometric diversity by widening the transducer pool area.
% - Increased the number of acquisitions to provide more data for reconstruction.
% - Optimized the imaging grid resolution to match the physical system.
% - Provided options for tuning the ADMM regularization parameter.
% =========================================================================
clear; clc; close all;

%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== UNIFIED PMUT ULTRASOUND SIMULATION ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('unified_simulation_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- Core Physical and Simulation Parameters ---
params = struct();
params.c = 343;                    % Speed of sound (m/s)
params.fs = 1e6;                   % Sampling frequency (Hz)
params.pmut_width_m = 0.020;       % pMUT width (m)
params.excitation_amplitude = 1e15;% High amplitude for strong signal
params.target_SNR_db = 35;         % Target SNR for adding noise to measurements

% --- H-Matrix Generation Parameters ---
% IMPROVEMENT: More acquisitions provide more information for reconstruction.
params.num_acquisitions = 50;      % (Original: 25)
params.num_active_tx = 5;          % Number of active transmitters per acquisition
params.max_delay_rand_us = 500;    % Max random delay for beamforming
params.apodization_mode = 'uniform'; % 'uniform' or 'random'
% IMPROVEMENT: A wider pool creates more diverse angles of incidence.
params.tx_pool_width_m = 0.300;    % (Original: 0.200)

% --- Realistic pMUT Waveform Parameters ---
params.pmut_resonance_freq = 57700; % 57.7 kHz average resonance
params.pmut_bandwidth = 2520;      % 2.52 kHz standard deviation of resonance
params.impulse_duration_us = 10;   % Short impulse excitation (10 us)
params.frequency_offset_hz = 0;    % Optional global offset from resonance

% --- Imaging Grid and Target Scene Parameters ---
params.grid_width_m = 0.150;       % Width of the imaging area
params.target_distance_m = 0.150;  % Center distance to imaging area
params.grid_depth_range_m = 0.100; % Depth of the imaging area
% IMPROVEMENT: A slightly coarser grid can prevent smearing and improve PSNR.
params.grid_step_m = 0.0075;       % (Original: 0.005)
params.target_pattern = '3x3_grid';% '3x3_grid', '2x2_grid', 'line_5', 'cross_5'
params.target_size_mm = 4;         % Diameter of targets in mm
params.grid_spacing_mm = 20;       % Spacing between targets in mm

% --- ADMM Reconstruction Parameters ---
params.numItersADMM = 50;
params.rho_admm = 6.73;
% IMPROVEMENT: Tuning lambda is crucial for reconstruction quality.
params.lambda_tv_reg = 1.2;      % (Original: 0.7438) Try values between 0.5 and 2.0
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;

%% ===== 2. RUN FULL SIMULATION =====

% --- Step 1: Generate the Sensing Matrix (H) ---
fprintf('\n--- Step 1: Generating H-Matrix ---\n');
tic;
[H, imaging_grid] = generate_h_matrix(params);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

% --- Step 2: Create the Target Scene (v) ---
fprintf('\n--- Step 2: Creating Target Scene ---\n');
tic;
scene_matrix = create_scene(imaging_grid, params);
v_true_vector = scene_matrix(:);
fprintf('Scene creation complete. Time: %.2f seconds.\n', toc);
visualize_scene(scene_matrix, imaging_grid, output_folder);

% --- Step 3: Run Forward Simulation (b = H*v + n) ---
fprintf('\n--- Step 3: Running Forward Simulation ---\n');
tic;
b_vector = run_forward_simulation(H, v_true_vector, params);
fprintf('Forward simulation complete. Time: %.2f seconds.\n', toc);

% --- Step 4: Reconstruct Image using ADMM ---
fprintf('\n--- Step 4: Reconstructing Image via ADMM ---\n');
tic;
reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);

% --- Step 5: Analyze and Save Results ---
fprintf('\n--- Step 5: Analyzing and Saving Results ---\n');
analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, output_folder);

fprintf('\n=== SIMULATION COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

% =========================================================================
% H-MATRIX GENERATION
% =========================================================================
function [H, imaging_grid] = generate_h_matrix(config)
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);
    
    vgrid_N = 10;
    vgrid_total_elements = vgrid_N * vgrid_N;
    vgrid_pitch = config.tx_pool_width_m / (vgrid_N - 1);
    vgrid_kerf = vgrid_pitch - config.pmut_width_m;
    if vgrid_kerf < 0.0001
        error('pMUT width is too large for the virtual grid.');
    end
    
    TxAperture = xdc_2d_array(vgrid_N, vgrid_N, config.pmut_width_m, config.pmut_width_m, vgrid_kerf, vgrid_kerf, ones(vgrid_N), 1, 1, [0 0 0]);
    RxAperture = xdc_2d_array(1, 1, config.pmut_width_m, config.pmut_width_m, 0, 0, ones(1,1), 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));

    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    
    imaging_grid = struct('X_mesh', X_mesh, 'Z_mesh', Z_mesh, 'points', grid_points);

    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        
        individual_resonances = config.pmut_resonance_freq + config.frequency_offset_hz + (rand(1, num_active_tx) - 0.5) * config.pmut_bandwidth;
        impulse_duration_samples = round(config.impulse_duration_us * fs / 1e6);
        max_len = impulse_duration_samples;
        
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            t = 0:1/fs:(impulse_duration_samples-1)/fs;
            random_phase = 2 * pi * rand();
            impulse_signal = sin(2 * pi * individual_resonances(k) * t + random_phase);
            window = tukeywin(length(t), 0.25)';
            tx_signal = impulse_signal .* window;
            composite_waveform(1:length(tx_signal)) = composite_waveform(1:length(tx_signal)) + tx_signal;
        end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        
        full_apod_vector = zeros(1, vgrid_total_elements);
        full_excitation_vector = zeros(1, vgrid_total_elements);
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        full_apod_vector(active_indices) = apod_weights;
        full_excitation_vector(active_indices) = 1;
        xdc_apodization(TxAperture, 0, full_apod_vector);
        xdc_excitation(TxAperture, full_excitation_vector);
        
        full_delay_vector = zeros(1, vgrid_total_elements);
        delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        full_delay_vector(active_indices) = delays_us * 1e-6;
        xdc_focus_times(TxAperture, 0, full_delay_vector);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    disp('  Assembling H-matrix using interpolation...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(config.num_acquisitions, N_pixels); return;
    end
    
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0, K_global_per_acq = 1; end
    
    total_rows = K_global_per_acq * config.num_acquisitions;
    estimated_nnz = round(sum(all_K_values) * N_pixels * 0.1);
    H = spalloc(total_rows, N_pixels, estimated_nnz);
    
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
    
    xdc_free(TxAperture);
    xdc_free(RxAperture);
    field_end();
end

% =========================================================================
% SCENE AND SIMULATION
% =========================================================================
function scene_matrix = create_scene(imaging_grid, params)
    scene_matrix = zeros(size(imaging_grid.X_mesh));
    X_mm = imaging_grid.X_mesh * 1000;
    Z_mm = imaging_grid.Z_mesh * 1000;
    
    grid_spacing_mm = params.grid_spacing_mm;
    target_size_mm = params.target_size_mm;
    grid_offset_x_mm = 0;
    grid_offset_z_mm = params.target_distance_m * 1000;
    
    switch params.target_pattern
        case '3x3_grid'
            positions = [];
            for row = -1:1, for col = -1:1
                positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm, 1.0];
            end, end
        case '2x2_grid'
            positions = [];
            for row = -0.5:1:0.5, for col = -0.5:1:0.5
                positions = [positions; col*grid_spacing_mm, row*grid_spacing_mm, 1.0];
            end, end
        case 'line_5'
            positions = ((-2:2)' * grid_spacing_mm);
            positions = [positions, zeros(5,1), ones(5,1)];
        case 'cross_5'
            positions = [0,0,1; -1,0,1; 1,0,1; 0,-1,1; 0,1,1] * grid_spacing_mm;
            positions(:,3) = 1;
    end
    
    target_radius_pixels = round(target_size_mm / (2 * (params.grid_step_m*1000)));
    if target_radius_pixels < 1, target_radius_pixels = 1; end
    
    for i = 1:size(positions, 1)
        x_pos_mm = positions(i, 1) + grid_offset_x_mm;
        z_pos_mm = positions(i, 2) + grid_offset_z_mm;
        
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iz_scene] = min(abs(Z_mm(:,1) - z_pos_mm));
        
        for dz = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                ix_target = ix_scene + dx;
                iz_target = iz_scene + dz;
                if ix_target > 0 && ix_target <= size(X_mm, 2) && iz_target > 0 && iz_target <= size(X_mm, 1)
                    scene_matrix(iz_target, ix_target) = positions(i, 3);
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

% =========================================================================
% ADMM RECONSTRUCTION
% =========================================================================
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
    for k_admm = 1:config.admm_max_iter
        x_prev = x_admm_img_iter;
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = admm_iteration(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            config.rho_admm, config.lambda_tv_reg, Hfun_pcg_admm, config);
        
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol
                fprintf('    ADMM converged at iteration %d.\n', k_admm);
                break;
            end
        end
    end
    reconstructed_image = x_admm_img_iter;
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [Dx, Dy] = difference_operators(imageSize)
   rows = imageSize(1); cols = imageSize(2);
   N_img_pixels = rows * cols;
   Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
   Dx( (cols-1)*rows+1 : cols*rows , :) = 0;
   Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
   Dy( rows:rows:N_img_pixels , :) = 0;
end

function [x_new, z_new, u_new] = admm_iteration(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, Hfun_pcg_admm, config)
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm = max(v_norm, eps);
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    
    u_new = u_old + opDx_tv(x_new) - z_new;
end

% =========================================================================
% ANALYSIS AND VISUALIZATION
% =========================================================================
function visualize_scene(scene_matrix, imaging_grid, output_folder)
    figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
    imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Z_mesh(:,1)*1000, scene_matrix);
    colormap(gray); colorbar; axis image;
    title('Ground Truth Target Scene');
    xlabel('X Position (mm)'); ylabel('Z Position (mm)');
    saveas(gcf, fullfile(output_folder, 'ground_truth_scene.png'));
    close(gcf);
end

function analyze_and_plot_results(reconstructed_image, scene_matrix, H, params, output_folder)
    % --- Calculate Metrics ---
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(abs(reconstructed_image(:)));
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    Hn = Hn ./ vecnorm(Hn, 2, 1);
    coherence_matrix = abs(Hn' * Hn);
    coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
    max_coherence = full(max(coherence_matrix(:)));
    
    fprintf('  PSNR: %.2f dB\n', psnr);
    fprintf('  Max Coherence: %.4f\n', max_coherence);
    
    % --- Create Plots ---
    figure('Visible', 'off', 'Position', [100, 100, 1200, 500]);
    
    % Ground Truth
    subplot(1, 2, 1);
    imagesc(scene_matrix); colormap(gray); axis image; colorbar;
    title('Ground Truth');
    xlabel('X Index'); ylabel('Z Index');
    
    % Reconstruction
    subplot(1, 2, 2);
    imagesc(reconstructed_image); colormap(gray); axis image; colorbar;
    title(sprintf('Reconstruction (PSNR: %.2f dB)', psnr));
    xlabel('X Index'); ylabel('Z Index');
    
    sgtitle('Simulation Results');
    saveas(gcf, fullfile(output_folder, 'reconstruction_comparison.png'));
    close(gcf);
    
    % --- Save Metrics ---
    metrics = struct('psnr', psnr, 'max_coherence', max_coherence);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
end
