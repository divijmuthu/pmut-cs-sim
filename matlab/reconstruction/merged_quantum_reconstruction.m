%% ===== MERGED QUANTUM RECONSTRUCTION =====
% Merges successful sweep script H matrix generation with TheoryTest2 reconstruction
% and improved target setup for optimal performance

clear; clc; close all;

% Add Field II path
addpath('m_files');

%% ===== OUTPUT SETUP =====
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('merged_quantum_output', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('=== MERGED QUANTUM RECONSTRUCTION ===\n');
fprintf('Combining sweep script H generation + TheoryTest2 reconstruction + improved targets\n');
fprintf('Saving results to: %s\n\n', output_folder);

%% ===== PARAMETERS (FROM VERIFIED SWEEP SCRIPTS) =====
% Physical parameters (from V29 sweep script)
params = struct();
params.c = 343;                    % Speed of sound (m/s) - from V29
params.fs = 1e6;                   % Sampling frequency (Hz) - from V29
params.pmut_width_m = 0.020;       % pMUT width (m) - from V29
params.tx_pool_width_m = 0.200;    % Transmitter pool width (m) - from V29
params.grid_width_m = 0.150;       % Imaging grid width (m) - from V29
params.target_distance_m = 0.150;  % Target distance (m) - from V29
params.grid_depth_range_m = 0.100; % Grid depth range (m) - from V29
params.grid_step_m = 0.010;        % Grid step (m) - from V29
params.num_acquisitions = 20;      % Number of acquisitions - from V29
params.excitation_amplitude = 1e15; % Excitation amplitude - from V29

% REALISTIC pMUT PARAMETERS (from experimental data)
params.pmut_resonance_freq = 57700; % 57.7 kHz average resonance - from V29
params.pmut_bandwidth = 2520;      % 2.52 kHz standard deviation - from V29
params.impulse_duration_us = 10;   % Short impulse excitation (10 μs) - from V29

% Fixed parameters (best from V28_final)
params.num_active_tx = 5;          % 5 transmitters (best from V28_final)
params.max_delay_rand_us = 500;    % 500μs delays (best from V28_final)
params.apodization_mode = 'uniform'; % Uniform apodization (best from V28_final)
params.frequency_offset_hz = 0;    % No frequency offset

% ADMM parameters (from TheoryTest2 - proven to work)
params.numItersADMM = 50;           % Fixed iterations
params.rho_admm = 6.73;              % Optimized ADMM penalty
params.lambda_tv_reg = 0.7438;       % Optimized TV regularization
params.admm_tol = 1.2e-5;             % Optimized tolerance
params.admm_max_iter = 50;          % Fixed max iterations
params.pcg_max_iter = 30;           % Reduced PCG iterations
params.pcg_tol = 1e-8;              % Slightly relaxed PCG tolerance

fprintf('Parameters configured:\n');
fprintf('  - pMUT resonance frequency: %.1f kHz\n', params.pmut_resonance_freq/1e3);
fprintf('  - Speed of sound: %d m/s\n', params.c);
fprintf('  - Grid width: %.1f m, depth range: %.1f m\n', params.grid_width_m, params.grid_depth_range_m);
fprintf('  - Number of acquisitions: %d\n', params.num_acquisitions);

%% ===== TARGET CONFIGURATIONS =====
target_configs = struct();

% High Challenge Configuration
target_configs.high_challenge = struct();
target_configs.high_challenge.target_size_mm = 2;
target_configs.high_challenge.grid_spacing_mm = 15;
target_configs.high_challenge.lambda_tv_reg = 0.1;
target_configs.high_challenge.description = 'High Challenge: 2mm targets, 15mm spacing';

% Optimal Challenge Configuration
target_configs.optimal = struct();
target_configs.optimal.target_size_mm = 4;
target_configs.optimal.grid_spacing_mm = 18;
target_configs.optimal.lambda_tv_reg = 0.1;
target_configs.optimal.description = 'Optimal Challenge: 4mm targets, 18mm spacing';

% Realistic Challenge Configuration
target_configs.realistic = struct();
target_configs.realistic.target_size_mm = 5;
target_configs.realistic.grid_spacing_mm = 20;
target_configs.realistic.lambda_tv_reg = 0.1;
target_configs.realistic.description = 'Realistic Challenge: 5mm targets, 20mm spacing';

%% ===== MAIN PROCESSING LOOP =====
config_names = fieldnames(target_configs);
num_configs = length(config_names);

for config_idx = 1:num_configs
    config_name = config_names{config_idx};
    config = target_configs.(config_name);
    
    fprintf('\n=== PROCESSING CONFIGURATION: %s ===\n', upper(config_name));
    fprintf('Description: %s\n', config.description);
    
    % Step 1: Create target scene
    fprintf('Creating target scene...\n');
    [target_scene, target_positions, imaging_grid] = create_improved_target_scene(config, params);
    
    % Step 2: Generate H matrix using sweep script approach
    fprintf('Generating H matrix using sweep script approach...\n');
    [H_matrix, transducer_positions] = generate_sweep_h_matrix(params, imaging_grid);
    
    % Step 3: Simulate measurements
    fprintf('Simulating measurements...\n');
    [measurements, measurement_times] = simulate_measurements(target_scene, H_matrix, params);
    
    % Step 4: Perform ADMM reconstruction using TheoryTest2 approach
    fprintf('Performing ADMM reconstruction...\n');
    [reconstructed_image, reconstruction_time, reconstruction_info] = perform_admm_reconstruction(measurements, H_matrix, config, params, imaging_grid);
    
    % Step 5: Calculate metrics
    fprintf('Calculating reconstruction metrics...\n');
    metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image);
    
    % Step 6: Visualize results
    fprintf('Creating visualizations...\n');
    visualize_merged_results(target_scene, reconstructed_image, H_matrix, measurements, ...
        target_positions, transducer_positions, imaging_grid, config, metrics, ...
        reconstruction_info, output_folder, config_name);
    
    % Step 7: Save results
    fprintf('Saving results...\n');
    save(fullfile(output_folder, sprintf('%s_merged_results.mat', config_name)), ...
        'target_scene', 'reconstructed_image', 'H_matrix', 'measurements', ...
        'target_positions', 'transducer_positions', 'imaging_grid', 'config', ...
        'metrics', 'reconstruction_time', 'reconstruction_info');
    
    fprintf('Configuration %s completed successfully!\n', config_name);
    fprintf('  - PSNR: %.2f dB\n', metrics.psnr);
    fprintf('  - Correlation: %.4f\n', metrics.correlation);
    fprintf('  - Reconstruction time: %.2f seconds\n', reconstruction_time);
end

%% ===== COMPARATIVE ANALYSIS =====
fprintf('\n=== COMPARATIVE ANALYSIS ===\n');
perform_comparative_analysis(output_folder, config_names);

fprintf('\n=== MERGED QUANTUM RECONSTRUCTION COMPLETE ===\n');
fprintf('All results saved to: %s\n', output_folder);

%% ===== HELPER FUNCTIONS =====

function [target_scene, target_positions, imaging_grid] = create_improved_target_scene(config, params)
    % Create improved target scene using V29 approach with high contrast
    
    % Calculate imaging grid (using V29 approach)
    x_coords_img = -params.grid_width_m/2 : params.grid_step_m : params.grid_width_m/2;
    z_coords_img = (params.target_distance_m - params.grid_depth_range_m/2) : params.grid_step_m : (params.target_distance_m + params.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    
    % Initialize target scene with background
    target_scene = ones(size(X_mesh)) * 0.1;  % Background intensity
    target_positions = [];
    
    % Create 3x3 grid of targets
    grid_center_x = 0;
    grid_center_z = params.target_distance_m;  % Use target distance from V29
    
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
            target_radius_pixels = max(1, round(target_size / params.grid_step_m));
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

function [H_matrix, transducer_positions] = generate_sweep_h_matrix(params, imaging_grid)
    % Generate H matrix using verified sweep script approach with improved coherence
    
    % Initialize Field II
    fprintf('  Initializing Field II...\n');
    field_init(-1);
    
    % Setup Field II transducers (using 2D array approach like V29)
    fs = params.fs;
    c = params.c;
    vgrid_N = 10;
    vgrid_total_elements = vgrid_N * vgrid_N;
    vgrid_pitch = params.tx_pool_width_m / (vgrid_N - 1);
    vgrid_kerf = vgrid_pitch - params.pmut_width_m;
    if vgrid_kerf < 0.0001
        error('pMUT width is too large for the virtual grid.');
    end
    
    fprintf('  Creating transducer arrays...\n');
    TxAperture = xdc_2d_array(vgrid_N, vgrid_N, params.pmut_width_m, params.pmut_width_m, vgrid_kerf, vgrid_kerf, ones(vgrid_N), 1, 1, [0 0 0]);
    RxAperture = xdc_2d_array(1, 1, params.pmut_width_m, params.pmut_width_m, 0, 0, ones(1,1), 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));

    % Create imaging grid points
    fprintf('  Creating imaging grid...\n');
    x_coords_img = -params.grid_width_m/2 : params.grid_step_m : params.grid_width_m/2;
    z_coords_img = (params.target_distance_m - params.grid_depth_range_m/2) : params.grid_step_m : (params.target_distance_m + params.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    
    fprintf('  Grid size: %d pixels\n', N_pixels);

    % Initialize H matrix with better coherence properties
    H_matrix = zeros(params.num_acquisitions, N_pixels);
    
    fprintf('  Starting %d acquisitions...\n', params.num_acquisitions);
    total_tic = tic;
    
    for r_acq = 1:params.num_acquisitions
        acq_tic = tic;
        
        % Use fixed number of transmitters (best from V28_final)
        num_active_tx = params.num_active_tx;
        
        fprintf('    Acquisition %d/%d: Using %d transmitters...', r_acq, params.num_acquisitions, num_active_tx);
        
        % Generate REALISTIC pMUT excitation (impulse at resonant frequency)
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        
        % REALISTIC: Each pMUT has slightly different resonant frequency (from experimental data)
        individual_resonances = params.pmut_resonance_freq + params.frequency_offset_hz + ...
            (rand(1, num_active_tx) - 0.5) * params.pmut_bandwidth;
        
        % REALISTIC: Generate impulse excitation for each pMUT
        impulse_duration_samples = round(params.impulse_duration_us * fs / 1e6);
        max_len = impulse_duration_samples;
        
        % Setup apodization
        apod_weights = ones(1, num_active_tx);
        if strcmp(params.apodization_mode, 'random')
            apod_weights = rand(1, num_active_tx);
        end
        
        % REALISTIC: Generate impulse excitation at resonant frequency
        excitation_amps = (0.5 + rand(1, num_active_tx)) * params.excitation_amplitude;
        composite_waveform = zeros(1, max_len);
        
        for k = 1:num_active_tx
            % REALISTIC: Short impulse at resonant frequency
            t = 0:1/fs:(impulse_duration_samples-1)/fs;
            random_phase = 2 * pi * rand();
            
            % Impulse excitation at pMUT's resonant frequency
            impulse_signal = sin(2 * pi * individual_resonances(k) * t + random_phase);
            
            % Apply short window (realistic impulse)
            window = ones(1, length(t));
            window(1:round(length(t)*0.1)) = linspace(0, 1, round(length(t)*0.1)); % Rise
            window(end-round(length(t)*0.1)+1:end) = linspace(1, 0, round(length(t)*0.1)); % Fall
            
            tx_signal = impulse_signal .* window * excitation_amps(k);
            composite_waveform(1:length(tx_signal)) = composite_waveform(1:length(tx_signal)) + tx_signal;
        end
        
        % Apply to Field II
        xdc_impulse(TxAperture, composite_waveform);
        
        % Setup apodization and excitation
        full_apod_vector = zeros(1, vgrid_total_elements);
        full_excitation_vector = zeros(1, vgrid_total_elements);
        full_apod_vector(active_indices) = apod_weights;
        full_excitation_vector(active_indices) = 1;
        xdc_apodization(TxAperture, 0, full_apod_vector);
        xdc_excitation(TxAperture, full_excitation_vector);
        
        % Setup delays
        full_delay_vector = zeros(1, vgrid_total_elements);
        if strcmp(params.apodization_mode, 'uniform')
            delays_us = zeros(1, num_active_tx);
        else
            delays_us = rand(1, num_active_tx) * params.max_delay_rand_us;
        end
        full_delay_vector(active_indices) = delays_us * 1e-6;
        xdc_focus_times(TxAperture, 0, full_delay_vector);
        
        % Calculate impulse response using Field II
        [h_r, start_time] = calc_hhp(TxAperture, RxAperture, grid_points);
        
        % Extract H matrix row with improved coherence properties
        if ~isempty(h_r) && size(h_r, 1) > 0
            % Use weighted sum over time samples to improve coherence
            time_weights = exp(-(1:size(h_r, 1)) / (size(h_r, 1) / 3)); % Exponential decay
            time_weights = time_weights / sum(time_weights); % Normalize
            H_matrix(r_acq, :) = sum(h_r .* time_weights', 1);
        end
        
        acq_time = toc(acq_tic);
        fprintf(' completed in %.2f seconds\n', acq_time);
    end
    
    total_time = toc(total_tic);
    fprintf('  All acquisitions completed in %.2f seconds\n', total_time);
    
    % Cleanup Field II
    xdc_free(TxAperture);
    xdc_free(RxAperture);
    field_end();
    
    % Calculate transducer positions for visualization
    transducer_width = vgrid_N * vgrid_pitch;
    transducer_positions = linspace(-transducer_width/2, transducer_width/2, vgrid_N);
    
    % Normalize H matrix to improve coherence
    column_norms = sqrt(sum(H_matrix.^2, 1));
    valid_columns = column_norms > 1e-10;
    H_matrix = H_matrix(:, valid_columns);
    
    % Apply coherence improvement techniques
    H_matrix = improve_coherence(H_matrix);
    
    fprintf('  H matrix generation completed\n');
    fprintf('  - Matrix dimensions: %d x %d\n', size(H_matrix, 1), size(H_matrix, 2));
end

function H_improved = improve_coherence(H)
    % Apply techniques to improve coherence while maintaining matrix properties
    
    % 1. Normalize columns
    H_norm = H ./ sqrt(sum(H.^2, 1));
    
    % 2. Apply slight randomization to break perfect correlations
    noise_level = 0.01; % 1% noise
    H_noisy = H_norm + noise_level * randn(size(H_norm));
    
    % 3. Re-normalize
    H_improved = H_noisy ./ sqrt(sum(H_noisy.^2, 1));
    
    % 4. Ensure no all-zero columns
    column_norms = sqrt(sum(H_improved.^2, 1));
    valid_columns = column_norms > 1e-10;
    H_improved = H_improved(:, valid_columns);
    
    fprintf('  Applied coherence improvement techniques\n');
end

function H = assemble_time_domain_h_matrix(all_h_data, all_start_times, all_K_values, params, N_pixels)
    % Assemble H matrix using time-domain interpolation (from successful sweep scripts)
    
    % Find valid acquisitions
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(params.num_acquisitions, N_pixels);
        fprintf('  WARNING: No valid acquisitions found!\n');
        return;
    end
    
    % Calculate global time window
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / params.fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/params.fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0, K_global_per_acq = 1; end
    
    fprintf('  Global time window: %.3f to %.3f ms (%d samples per acquisition)\n', ...
        min_global_start_time*1000, max_global_end_time*1000, K_global_per_acq);
    
    % Pre-allocate sparse matrix
    total_rows = K_global_per_acq * params.num_acquisitions;
    estimated_nnz = round(sum(all_K_values) * N_pixels * 0.1);
    H = spalloc(total_rows, N_pixels, estimated_nnz);
    
    % Assemble H matrix with time-domain interpolation
    current_row_offset = 0;
    
    for r_acq = 1:params.num_acquisitions
        h_r = all_h_data{r_acq};
        start_time = all_start_times(r_acq);
        K_current = all_K_values(r_acq);
        
        if K_current == 0 || isempty(h_r)
            current_row_offset = current_row_offset + K_global_per_acq;
            continue;
        end
        
        % Calculate time axis for this acquisition
        t_current_acq_axis = start_time + (0:(K_current - 1)) / params.fs;
        
        % Interpolate to common time axis
        h_interpolated = interpolate_to_common_time(t_current_acq_axis, h_r, t_common_axis, K_global_per_acq, N_pixels);
        
        % Assign to H matrix
        row_indices = current_row_offset + (1:K_global_per_acq);
        if max(row_indices) <= size(H, 1)
            H(row_indices, :) = h_interpolated;
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    
    % Convert to sparse matrix
    H = sparse(H);
    
    fprintf('  Time-domain H matrix assembly completed\n');
    fprintf('  - Final matrix: %d x %d\n', size(H, 1), size(H, 2));
    fprintf('  - Sparsity: %.2f%%\n', 100 * (1 - nnz(H) / numel(H)));
end

function h_interpolated = interpolate_to_common_time(t_current_acq_axis, h_current, t_common_axis, K_global_per_acq, N_pixels)
    % Interpolate acquisition data to common time axis
    
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

function [measurements, measurement_times] = simulate_measurements(target_scene, H_matrix, params)
    % Simulate measurements
    
    % Flatten target scene
    target_scene_flat = target_scene(:);
    
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

function [reconstructed_image, reconstruction_time, reconstruction_info] = perform_admm_reconstruction(measurements, H_matrix, config, params, imaging_grid)
    % Perform ADMM reconstruction using TheoryTest2 approach
    
    tic;
    
    % Calculate imaging grid dimensions
    num_x = length(imaging_grid.x_grid);
    num_z = length(imaging_grid.z_grid);
    
    % Initialize reconstruction
    num_pixels = size(H_matrix, 2);
    
    % ADMM parameters (from TheoryTest2)
    lambda_tv = config.lambda_tv_reg;
    rho = params.rho_admm;
    max_iter = params.admm_max_iter;
    tol = params.admm_tol;
    
    % Initialize variables
    x = zeros(num_pixels, 1);  % Main variable
    z = zeros(num_pixels, 1);  % Auxiliary variable
    u = zeros(num_pixels, 1);  % Dual variable
    
    % Precompute matrices
    HtH = H_matrix' * H_matrix;
    Hty = H_matrix' * measurements;
    
    % ADMM iterations
    reconstruction_info.iterations = 0;
    reconstruction_info.converged = false;
    
    for iter = 1:max_iter
        % Update x (least squares step)
        x_old = x;
        x = (HtH + rho * eye(num_pixels)) \ (Hty + rho * (z - u));
        
        % Check for NaN or Inf values
        if any(isnan(x)) || any(isinf(x))
            fprintf('  Warning: NaN/Inf detected in x at iteration %d, using simple least squares\n', iter);
            x = H_matrix \ measurements;
            break;
        end
        
        % Update z (TV regularization step)
        z_old = z;
        z = soft_threshold_tv(x + u, lambda_tv / rho);
        
        % Update dual variable
        u = u + x - z;
        
        % Check convergence
        primal_residual = norm(x - z);
        dual_residual = rho * norm(z - z_old);
        
        if primal_residual < tol && dual_residual < tol
            reconstruction_info.converged = true;
            reconstruction_info.iterations = iter;
            break;
        end
    end
    
    % Reshape to image
    reconstructed_image = reshape(x, num_z, num_x);
    
    reconstruction_time = toc;
    
    % Store reconstruction info
    reconstruction_info.final_primal_residual = primal_residual;
    reconstruction_info.final_dual_residual = dual_residual;
    reconstruction_info.lambda_tv = lambda_tv;
    reconstruction_info.rho = rho;
    
    fprintf('  ADMM reconstruction completed in %.2f seconds\n', reconstruction_time);
    fprintf('  - Iterations: %d\n', reconstruction_info.iterations);
    fprintf('  - Converged: %s\n', mat2str(reconstruction_info.converged));
    fprintf('  - Final primal residual: %.2e\n', reconstruction_info.final_primal_residual);
end

function z = soft_threshold_tv(x, threshold)
    % Simplified soft thresholding for TV regularization
    
    % Apply simple soft thresholding to the vector
    z = soft_threshold(x, threshold);
    
    % Apply basic smoothing
    z = smooth(z, 3);
end

function y = soft_threshold(x, threshold)
    % Soft thresholding function
    y = sign(x) .* max(abs(x) - threshold, 0);
end

function y = smooth(x, window_size)
    % Simple smoothing function
    if length(x) < window_size
        y = x;
        return;
    end
    
    y = x;
    for i = 1:length(x)
        start_idx = max(1, i - floor(window_size/2));
        end_idx = min(length(x), i + floor(window_size/2));
        y(i) = mean(x(start_idx:end_idx));
    end
end

function metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image)
    % Calculate reconstruction quality metrics with meaningful comparison
    
    fprintf('  === METRICS DEBUGGING ===\n');
    fprintf('  Target scene: min=%.6f, max=%.6f, mean=%.6f, std=%.6f\n', ...
        min(target_scene(:)), max(target_scene(:)), mean(target_scene(:)), std(target_scene(:)));
    fprintf('  Reconstructed image: min=%.6f, max=%.6f, mean=%.6f, std=%.6f\n', ...
        min(reconstructed_image(:)), max(reconstructed_image(:)), mean(reconstructed_image(:)), std(reconstructed_image(:)));
    
    % Check for all-zero or constant images
    if max(target_scene(:)) == min(target_scene(:))
        fprintf('  WARNING: Target scene is constant!\n');
    end
    if max(reconstructed_image(:)) == min(reconstructed_image(:))
        fprintf('  WARNING: Reconstructed image is constant!\n');
    end
    
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
    
    fprintf('  After normalization:\n');
    fprintf('  Target norm: min=%.6f, max=%.6f, mean=%.6f\n', ...
        min(target_norm(:)), max(target_norm(:)), mean(target_norm(:)));
    fprintf('  Recon norm: min=%.6f, max=%.6f, mean=%.6f\n', ...
        min(recon_norm(:)), max(recon_norm(:)), mean(recon_norm(:)));
    
    % Calculate PSNR using the normalized images
    mse = mean((target_norm(:) - recon_norm(:)).^2);
    fprintf('  MSE: %.10f\n', mse);
    
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
    
    fprintf('  === FINAL METRICS ===\n');
    fprintf('  PSNR: %.2f dB, Correlation: %.4f, SSIM: %.4f, Rel Error: %.4f\n', ...
        metrics.psnr, metrics.correlation, metrics.ssim, metrics.relative_error);
    fprintf('  === END DEBUGGING ===\n');
end

function visualize_merged_results(target_scene, reconstructed_image, H_matrix, measurements, ...
    target_positions, transducer_positions, imaging_grid, config, metrics, ...
    reconstruction_info, output_folder, config_name)
    % Create comprehensive visualization of merged reconstruction results
    
    figure('Position', [100, 100, 1800, 1200]);
    
    % 1. Original target scene
    subplot(4, 5, 1);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, target_scene);
    colorbar;
    xlabel('X (m)'); ylabel('Z (m)');
    title('Original Target Scene');
    axis equal tight;
    
    % 2. Raw reconstructed image
    subplot(4, 5, 2);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, reconstructed_image);
    colorbar;
    xlabel('X (m)'); ylabel('Z (m)');
    title('Raw Reconstructed Image');
    axis equal tight;
    
    % 3. Normalized reconstructed image
    subplot(4, 5, 3);
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, recon_norm);
    colorbar;
    xlabel('X (m)'); ylabel('Z (m)');
    title('Normalized Reconstructed Image');
    axis equal tight;
    caxis([0 1]);
    
    % 4. Difference image
    subplot(4, 5, 4);
    diff_image = abs(target_scene - reconstructed_image);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, diff_image);
    colorbar;
    xlabel('X (m)'); ylabel('Z (m)');
    title('Absolute Difference');
    axis equal tight;
    
    % 5. Target positions overlay
    subplot(4, 5, 5);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, target_scene);
    hold on;
    scatter(target_positions(:,1), target_positions(:,2), 100, 'r', 'filled');
    for i = 1:size(target_positions, 1)
        text(target_positions(i,1), target_positions(i,2), sprintf('T%d', i), ...
            'Color', 'white', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
    xlabel('X (m)'); ylabel('Z (m)');
    title('Target Positions');
    axis equal tight;
    
    % 6. H matrix visualization
    subplot(4, 5, 6);
    imagesc(H_matrix);
    colorbar;
    xlabel('Pixel Index'); ylabel('Acquisition Index');
    title('H Matrix');
    
    % 7. Measurements
    subplot(4, 5, 7);
    plot(measurements);
    xlabel('Measurement Index'); ylabel('Amplitude');
    title('Measurements');
    grid on;
    
    % 8. Transducer positions
    subplot(4, 5, 8);
    scatter(transducer_positions, zeros(size(transducer_positions)), 50, 'b', 'filled');
    xlabel('X (m)'); ylabel('Y (m)');
    title('Transducer Array');
    axis equal; grid on;
    
    % 9. Metrics summary
    subplot(4, 5, 9);
    metrics_text = sprintf(['Reconstruction Metrics:\n\n' ...
                           'PSNR: %.2f dB\n' ...
                           'Correlation: %.4f\n' ...
                           'SSIM: %.4f\n' ...
                           'Relative Error: %.4f\n\n' ...
                           'Configuration:\n' ...
                           'Target Size: %.1f mm\n' ...
                           'Grid Spacing: %.1f mm'], ...
                           metrics.psnr, metrics.correlation, metrics.ssim, metrics.relative_error, ...
                           config.target_size_mm, config.grid_spacing_mm);
    text(0.1, 0.5, metrics_text, 'FontSize', 10, 'VerticalAlignment', 'middle');
    axis off;
    
    % 10. Reconstruction info
    subplot(4, 5, 10);
    recon_text = sprintf(['Reconstruction Info:\n\n' ...
                         'Iterations: %d\n' ...
                         'Converged: %s\n' ...
                         'Primal Residual: %.2e\n' ...
                         'Dual Residual: %.2e\n' ...
                         'Lambda TV: %.4f\n' ...
                         'Rho: %.4f'], ...
                         reconstruction_info.iterations, mat2str(reconstruction_info.converged), ...
                         reconstruction_info.final_primal_residual, reconstruction_info.final_dual_residual, ...
                         reconstruction_info.lambda_tv, reconstruction_info.rho);
    text(0.1, 0.5, recon_text, 'FontSize', 10, 'VerticalAlignment', 'middle');
    axis off;
    
    sgtitle(sprintf('Merged Quantum Reconstruction: %s Configuration', upper(config_name)), 'FontSize', 16);
    
    % Save figure
    saveas(gcf, fullfile(output_folder, sprintf('%s_merged_results.png', config_name)));
    fprintf('  Visualization saved: %s_merged_results.png\n', config_name);
end

function perform_comparative_analysis(output_folder, config_names)
    % Perform comparative analysis of all configurations
    
    fprintf('Loading results for comparative analysis...\n');
    
    comparative_results = struct();
    for config_idx = 1:length(config_names)
        config_name = config_names{config_idx};
        results_file = fullfile(output_folder, sprintf('%s_merged_results.mat', config_name));
        
        if exist(results_file, 'file')
            load(results_file);
            comparative_results.(config_name) = struct();
            comparative_results.(config_name).metrics = metrics;
            comparative_results.(config_name).reconstruction_time = reconstruction_time;
            comparative_results.(config_name).config = config;
        else
            fprintf('  Warning: Could not find results for %s\n', config_name);
        end
    end
    
    % Create comparative visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % PSNR comparison
    subplot(2, 3, 1);
    psnr_values = [];
    config_labels = {};
    for config_idx = 1:length(config_names)
        config_name = config_names{config_idx};
        if isfield(comparative_results, config_name)
            psnr_values = [psnr_values, comparative_results.(config_name).metrics.psnr];
            config_labels{end+1} = config_name;
        end
    end
    bar(psnr_values);
    set(gca, 'XTickLabel', config_labels);
    title('PSNR Comparison');
    ylabel('PSNR (dB)');
    grid on;
    
    % Correlation comparison
    subplot(2, 3, 2);
    corr_values = [];
    for config_idx = 1:length(config_names)
        config_name = config_names{config_idx};
        if isfield(comparative_results, config_name)
            corr_values = [corr_values, comparative_results.(config_name).metrics.correlation];
        end
    end
    bar(corr_values);
    set(gca, 'XTickLabel', config_labels);
    title('Correlation Comparison');
    ylabel('Correlation');
    grid on;
    
    % Reconstruction time comparison
    subplot(2, 3, 3);
    time_values = [];
    for config_idx = 1:length(config_names)
        config_name = config_names{config_idx};
        if isfield(comparative_results, config_name)
            time_values = [time_values, comparative_results.(config_name).reconstruction_time];
        end
    end
    bar(time_values);
    set(gca, 'XTickLabel', config_labels);
    title('Reconstruction Time Comparison');
    ylabel('Time (seconds)');
    grid on;
    
    % Summary table
    subplot(2, 3, [4, 5, 6]);
    summary_data = {};
    for config_idx = 1:length(config_names)
        config_name = config_names{config_idx};
        if isfield(comparative_results, config_name)
            summary_data{end+1, 1} = config_name;
            summary_data{end, 2} = sprintf('%.2f dB', comparative_results.(config_name).metrics.psnr);
            summary_data{end, 3} = sprintf('%.4f', comparative_results.(config_name).metrics.correlation);
            summary_data{end, 4} = sprintf('%.2f s', comparative_results.(config_name).reconstruction_time);
        end
    end
    
    if ~isempty(summary_data)
        column_names = {'Configuration', 'PSNR', 'Correlation', 'Time'};
        uitable('Data', summary_data, 'ColumnName', column_names, 'Position', [0.1, 0.1, 0.8, 0.8]);
    end
    
    sgtitle('Merged Quantum Reconstruction: Comparative Analysis', 'FontSize', 16);
    
    % Save comparative analysis
    saveas(gcf, fullfile(output_folder, 'comparative_analysis.png'));
    fprintf('  Comparative analysis saved: comparative_analysis.png\n');
    
    % Find best configuration
    if ~isempty(psnr_values)
        [best_psnr, best_idx] = max(psnr_values);
        fprintf('Best PSNR: %.2f dB (%s configuration)\n', best_psnr, config_labels{best_idx});
    end
end 