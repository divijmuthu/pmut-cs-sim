%% visualize_best_configuration.m - Visualize Best Configuration from Parameter Sweep
% Generate visualizations for the best performing configuration

clear; clc; close all;

%% ===== BEST CONFIGURATION (from parameter sweep results) =====
% Best configuration from our analysis:
% PSNR: 52.40 dB, Grid Step: 7.5 mm, Target Size: 3.0 mm, Grid Spacing: 20.0 mm

fprintf('=== VISUALIZING BEST CONFIGURATION ===\n');
fprintf('Best PSNR: 52.40 dB\n');
fprintf('Grid Step: 7.5 mm\n');
fprintf('Target Size: 3.0 mm\n');
fprintf('Grid Spacing: 20.0 mm\n\n');

%% ===== CONFIGURATION SETUP =====
config = struct();

% Physical parameters
config.c = 343;                    % Speed of sound (m/s)
config.fs = 1e6;                   % Sampling frequency (Hz)
config.pmut_width_m = 0.020;       % pMUT width (m)
config.tx_pool_width_m = 0.200;    % Transmitter pool width (m)
config.grid_width_m = 0.150;       % Imaging grid width (m)
config.target_distance_m = 0.150;  % Target distance (m)
config.grid_depth_range_m = 0.100; % Grid depth range (m)

% BEST PARAMETERS FROM SWEEP
config.grid_step_m = 0.0075;       % 7.5 mm grid step (BEST)
config.num_acquisitions = 20;      % Number of acquisitions
config.excitation_amplitude = 1e15; % Excitation amplitude

% Realistic pMUT parameters
config.pmut_resonance_freq = 57700; % 57.7 kHz average resonance
config.pmut_bandwidth = 2520;      % 2.52 kHz standard deviation
config.impulse_duration_us = 10;   % Short impulse excitation (10 μs)

% ADMM parameters
config.admm_max_iter = 50;
config.pcg_tol = 1e-8;
config.pcg_max_iter = 100;
config.rho_admm = 1.0;
config.lambda_tv_reg = 0.7438;
config.target_SNR_db = 35;

% Scene parameters (BEST FROM SWEEP)
config.target_size_mm = 3.0;       % 3.0 mm targets (BEST)
config.grid_spacing_mm = 20.0;     % 20.0 mm spacing (BEST)
config.target_pattern = '3x3_grid'; % 3x3 grid pattern

%% ===== OUTPUT SETUP =====
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('best_config_visualization', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('Saving visualizations to: %s\n\n', output_folder);

%% ===== GENERATE H MATRIX =====
fprintf('Generating H matrix with best configuration...\n');
tic;
[H, coherence_matrix] = generate_h_matrix_realistic_pmut(config);
h_generation_time = toc;
fprintf('H matrix generation completed in %.2f seconds\n', h_generation_time);

%% ===== CREATE IMAGING GRID =====
fprintf('Creating imaging grid...\n');
x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);

fprintf('Grid size: %d x %d = %d pixels\n', size(X_mesh, 1), size(X_mesh, 2), numel(X_mesh));

%% ===== CREATE SCENE =====
fprintf('Creating scene with best parameters...\n');
scene_matrix = create_scene_with_pattern(X_mesh, Z_mesh, config, config);
v_true_vector = scene_matrix(:);

% Plot true scene
plot_true_scene(scene_matrix, X_mesh, Z_mesh, config, output_folder);

%% ===== FORWARD SIMULATION =====
fprintf('Running forward simulation...\n');
tic;
Hv_signal = H * v_true_vector;
signal_power = mean(Hv_signal(:).^2);
target_SNR_linear = 10^(config.target_SNR_db / 10);
noise_variance = signal_power / target_SNR_linear;
noise_sigma = sqrt(noise_variance);
noise = noise_sigma * randn(size(Hv_signal));
b_vector = Hv_signal + noise;
forward_time = toc;

fprintf('Forward simulation completed in %.2f seconds\n', forward_time);
fprintf('Signal power: %.2e, Noise power: %.2e, SNR: %.1f dB\n', ...
    signal_power, noise_variance, 10*log10(signal_power/noise_variance));

% Plot measurement signals
plot_measurement_signals(b_vector, config, output_folder);

%% ===== ADMM RECONSTRUCTION =====
fprintf('Running ADMM reconstruction...\n');
tic;
reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, output_folder);
reconstruction_time = toc;

fprintf('Reconstruction completed in %.2f seconds\n', reconstruction_time);

%% ===== CALCULATE METRICS =====
fprintf('Calculating reconstruction metrics...\n');

scene_norm = scene_matrix / max(scene_matrix(:));
recon_norm = reconstructed_image / max(reconstructed_image(:));

MSE = mean((scene_norm(:) - recon_norm(:)).^2);
psnr = 10 * log10(1 / MSE);

scene_vec = scene_norm(:);
recon_vec = recon_norm(:);
correlation = (scene_vec - mean(scene_vec))' * (recon_vec - mean(recon_vec)) / ...
    (sqrt(sum((scene_vec - mean(scene_vec)).^2)) * sqrt(sum((recon_vec - mean(recon_vec)).^2)));

fprintf('\n=== RECONSTRUCTION METRICS ===\n');
fprintf('PSNR: %.2f dB\n', psnr);
fprintf('Correlation: %.4f\n', correlation);
fprintf('MSE: %.6f\n', MSE);

%% ===== PLOT RECONSTRUCTION COMPARISON =====
fprintf('Creating reconstruction comparison plot...\n');
plot_reconstruction_comparison(scene_matrix, reconstructed_image, X_mesh, Z_mesh, config, output_folder);

%% ===== SAVE RESULTS =====
fprintf('Saving results...\n');
save(fullfile(output_folder, 'best_config_results.mat'), 'config', 'psnr', 'correlation', 'MSE', ...
    'scene_matrix', 'reconstructed_image', 'H', 'b_vector', 'h_generation_time', 'forward_time', 'reconstruction_time');

fprintf('\n=== VISUALIZATION COMPLETE ===\n');
fprintf('Results saved to: %s\n', output_folder);
fprintf('Generated files:\n');
fprintf('  - true_scene.png (True target scene)\n');
fprintf('  - measurement_signals.png (Forward simulation signals)\n');
fprintf('  - admm_convergence.png (ADMM convergence)\n');
fprintf('  - reconstruction_comparison.png (True vs Reconstructed)\n');
fprintf('  - best_config_results.mat (Complete results)\n');

%% ===== HELPER FUNCTIONS =====

function [H, coherence_matrix] = generate_h_matrix_realistic_pmut(config)
    % Initialize Field II
    fprintf('  Initializing Field II...\n');
    field_init(-1);
    
    % Setup Field II transducers (using 2D array approach)
    fs = config.fs;
    c = config.c;
    vgrid_N = 10;
    vgrid_total_elements = vgrid_N * vgrid_N;
    vgrid_pitch = config.tx_pool_width_m / (vgrid_N - 1);
    vgrid_kerf = vgrid_pitch - config.pmut_width_m;
    if vgrid_kerf < 0.0001
        error('pMUT width is too large for the virtual grid.');
    end
    
    fprintf('  Creating transducer arrays...\n');
    TxAperture = xdc_2d_array(vgrid_N, vgrid_N, config.pmut_width_m, config.pmut_width_m, vgrid_kerf, vgrid_kerf, ones(vgrid_N), 1, 1, [0 0 0]);
    RxAperture = xdc_2d_array(1, 1, config.pmut_width_m, config.pmut_width_m, 0, 0, ones(1,1), 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));

    % Create imaging grid
    fprintf('  Creating imaging grid...\n');
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    fprintf('  Grid size: %d pixels\n', N_pixels);

    % Initialize acquisition storage
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    fprintf('  Starting %d acquisitions...\n', config.num_acquisitions);
    total_tic = tic;
    
    for r_acq = 1:config.num_acquisitions
        acq_tic = tic;
        
        % Use fixed number of transmitters (best from V28_final)
        num_active_tx = 5;
        
        fprintf('    Acquisition %d/%d: Using %d transmitters...', r_acq, config.num_acquisitions, num_active_tx);
        
        % Generate REALISTIC pMUT excitation (impulse at resonant frequency)
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        
        % REALISTIC: Each pMUT has slightly different resonant frequency (from experimental data)
        individual_resonances = config.pmut_resonance_freq + ...
            (rand(1, num_active_tx) - 0.5) * config.pmut_bandwidth;
        
        % REALISTIC: Generate impulse excitation for each pMUT
        impulse_duration_samples = round(config.impulse_duration_us * fs / 1e6);
        max_len = impulse_duration_samples;
        
        % Setup apodization (uniform)
        apod_weights = ones(1, num_active_tx);
        
        % REALISTIC: Generate impulse excitation at resonant frequency
        excitation_amps = (0.5 + rand(1, num_active_tx)) * config.excitation_amplitude;
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
        
        % Setup delays (uniform)
        full_delay_vector = zeros(1, vgrid_total_elements);
        delays_us = zeros(1, num_active_tx);
        full_delay_vector(active_indices) = delays_us * 1e-6;
        xdc_focus_times(TxAperture, 0, full_delay_vector);
        
        % Calculate impulse response
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        acq_time = toc(acq_tic);
        fprintf(' completed in %.2f seconds\n', acq_time);
    end
    
    total_time = toc(total_tic);
    fprintf('  All acquisitions completed in %.2f seconds\n', total_time);
    fprintf('  Assembling H-matrix using interpolation...\n');
    
    % Assemble H matrix using time-domain approach
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(config.num_acquisitions, N_pixels);
        coherence_matrix = [];
        fprintf('  WARNING: No valid acquisitions found!\n');
        return;
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
    
    fprintf('  H matrix size: %d x %d\n', total_rows, N_pixels);
    
    current_row_offset = 0;
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned = interpolate_h_response(t_current, all_h_data{r_acq}, t_common_axis);
            row_indices = current_row_offset + (1:K_global_per_acq);
            if max(row_indices) <= size(H, 1)
                H(row_indices, :) = h_aligned;
            end
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    
    % Cleanup Field II
    xdc_free(TxAperture);
    xdc_free(RxAperture);
    field_end();
    
    % H matrix statistics
    if nnz(H) > 0
        fprintf('  H matrix stats: min=%.3e, max=%.3e, nnz=%d, sum=%.3e\n', ...
            min(nonzeros(H)), max(nonzeros(H)), nnz(H), sum(nonzeros(H)));
    else
        fprintf('  H matrix stats: all zeros, nnz=0\n');
    end
    
    % Column analysis
    col_norms = vecnorm(H, 2, 1);
    num_nz_cols = sum(col_norms > 1e-20);
    fprintf('  Number of nonzero columns: %d (out of %d)\n', full(num_nz_cols), full(size(H,2)));
    
    % Compute coherence matrix for plotting
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    
    if isempty(Hn) || size(Hn, 2) < 2
        coherence_matrix = [];
    else
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0; % Remove diagonal
    end
end

function h_aligned = interpolate_h_response(t_current, h_current, t_common)
    % Interpolate H response to common time axis
    N_pixels = size(h_current, 2);
    h_aligned = zeros(length(t_common), N_pixels);
    if length(t_current) > 1
        for px_col = 1:N_pixels
            h_aligned(:, px_col) = interp1(t_current, h_current(:, px_col), t_common, 'linear', 0);
        end
    end
end

function scene_matrix = create_scene_with_pattern(X_mesh, Z_mesh, config, params)
    % Create scene based on pattern type
    
    scene_matrix = zeros(size(X_mesh));
    
    % Convert to mm for easier visualization
    X_mm = X_mesh * 1000;
    Z_mm = Z_mesh * 1000;
    
    % Extract scene parameters
    target_size_mm = params.target_size_mm;
    grid_spacing_mm = params.grid_spacing_mm;
    target_pattern = params.target_pattern;
    
    % Grid step in mm
    grid_step_mm = config.grid_step_m * 1000;
    target_radius_pixels = round(target_size_mm / (2 * grid_step_mm));
    
    % Ensure minimum size
    if target_radius_pixels < 1
        target_radius_pixels = 1;
    end
    
    % Create targets based on pattern
    if strcmp(target_pattern, '3x3_grid')
        % 3x3 grid of targets
        grid_offset_x_mm = 0;   % Center the grid
        grid_offset_z_mm = 150; % Position in middle of Z range
        
        target_positions = [];
        for row = 1:3
            for col = 1:3
                x_pos_mm = grid_offset_x_mm + (col - 2) * grid_spacing_mm;  % Center at 0
                z_pos_mm = grid_offset_z_mm + (row - 2) * grid_spacing_mm;  % Center at 150mm
                target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
            end
        end
        
    elseif strcmp(target_pattern, '2x2_grid')
        % 2x2 grid of targets
        grid_offset_x_mm = 0;
        grid_offset_z_mm = 150;
        
        target_positions = [];
        for row = 1:2
            for col = 1:2
                x_pos_mm = grid_offset_x_mm + (col - 1.5) * grid_spacing_mm;
                z_pos_mm = grid_offset_z_mm + (row - 1.5) * grid_spacing_mm;
                target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
            end
        end
        
    elseif strcmp(target_pattern, 'line_5')
        % Line of 5 targets
        grid_offset_x_mm = 0;
        grid_offset_z_mm = 150;
        
        target_positions = [];
        for i = 1:5
            x_pos_mm = grid_offset_x_mm + (i - 3) * grid_spacing_mm;
            z_pos_mm = grid_offset_z_mm;
            target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
        end
        
    elseif strcmp(target_pattern, 'cross_5')
        % Cross pattern with 5 targets
        grid_offset_x_mm = 0;
        grid_offset_z_mm = 150;
        
        target_positions = [
            grid_offset_x_mm, grid_offset_z_mm, 1.0;  % Center
            grid_offset_x_mm - grid_spacing_mm, grid_offset_z_mm, 1.0;  % Left
            grid_offset_x_mm + grid_spacing_mm, grid_offset_z_mm, 1.0;  % Right
            grid_offset_x_mm, grid_offset_z_mm - grid_spacing_mm, 1.0;  % Top
            grid_offset_x_mm, grid_offset_z_mm + grid_spacing_mm, 1.0;  % Bottom
        ];
        
    else
        error('Unknown target pattern: %s', target_pattern);
    end
    
    % Place targets
    for i = 1:size(target_positions, 1)
        x_pos_mm = target_positions(i, 1);
        z_pos_mm = target_positions(i, 2);
        amplitude = target_positions(i, 3);
        
        % Find closest grid points
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iz_scene] = min(abs(Z_mm(:,1) - z_pos_mm));
        
        % Create target
        for dz = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                % Check bounds
                if iz_scene + dz >= 1 && iz_scene + dz <= size(scene_matrix, 1) && ...
                   ix_scene + dx >= 1 && ix_scene + dx <= size(scene_matrix, 2)
                    scene_matrix(iz_scene + dz, ix_scene + dx) = amplitude;
                end
            end
        end
    end
    
    fprintf('  Scene created with %d targets:\n', size(target_positions, 1));
    fprintf('    Target size: %.1f mm\n', target_size_mm);
    fprintf('    Grid spacing: %.1f mm\n', grid_spacing_mm);
    fprintf('    Pattern: %s\n', target_pattern);
    
    % Print target positions
    for i = 1:size(target_positions, 1)
        fprintf('    Target %d: (%.1f, %.1f) mm, amplitude %.1f\n', ...
            i, target_positions(i, 1), target_positions(i, 2), target_positions(i, 3));
    end
end

function plot_true_scene(scene_matrix, X_mesh, Z_mesh, config, output_folder)
    % Plot the true scene
    
    figure('Position', [100, 100, 800, 600]);
    
    % Convert to mm for display
    X_mm = X_mesh * 1000;
    Z_mm = Z_mesh * 1000;
    
    % Create the plot
    imagesc(X_mm(1,:), Z_mm(:,1), scene_matrix);
    colormap(gray);
    colorbar;
    
    % Set axis properties
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    title(sprintf('True Scene: %s Pattern\nTarget Size: %.1f mm, Grid Spacing: %.1f mm', ...
        config.target_pattern, config.target_size_mm, config.grid_spacing_mm));
    
    % Add target markers
    hold on;
    [target_rows, target_cols] = find(scene_matrix > 0);
    if ~isempty(target_rows)
        target_x = X_mm(sub2ind(size(X_mm), target_rows, target_cols));
        target_z = Z_mm(sub2ind(size(Z_mm), target_rows, target_cols));
        scatter(target_x, target_z, 50, 'r', 'filled', 'MarkerEdgeColor', 'white');
    end
    hold off;
    
    % Save plot
    plot_filename = fullfile(output_folder, 'true_scene.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  True scene plot saved: %s\n', plot_filename);
end

function plot_measurement_signals(b_vector, config, output_folder)
    % Plot measurement signals
    
    figure('Position', [100, 100, 1000, 600]);
    
    % Reshape to show acquisitions
    signals_per_acq = length(b_vector) / config.num_acquisitions;
    b_matrix = reshape(b_vector, signals_per_acq, config.num_acquisitions);
    
    % Plot first few acquisitions
    num_to_plot = min(5, config.num_acquisitions);
    
    subplot(2, 1, 1);
    for i = 1:num_to_plot
        plot(b_matrix(:, i), 'LineWidth', 1.5);
        hold on;
    end
    xlabel('Time Sample');
    ylabel('Amplitude');
    title(sprintf('Measurement Signals (First %d Acquisitions)', num_to_plot));
    legend(arrayfun(@(x) sprintf('Acq %d', x), 1:num_to_plot, 'UniformOutput', false));
    grid on;
    
    subplot(2, 1, 2);
    % Show signal statistics
    signal_power = mean(b_vector.^2);
    noise_estimate = std(b_vector);
    snr_estimate = 20 * log10(signal_power / noise_estimate^2);
    
    histogram(b_vector, 50, 'Normalization', 'probability');
    xlabel('Signal Amplitude');
    ylabel('Probability');
    title(sprintf('Signal Distribution (SNR ≈ %.1f dB)', snr_estimate));
    grid on;
    
    % Save plot
    plot_filename = fullfile(output_folder, 'measurement_signals.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  Measurement signals plot saved: %s\n', plot_filename);
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, output_folder)
    % Run ADMM reconstruction with TV regularization
    
    fprintf('  Starting ADMM reconstruction...\n');
    
    % Setup ADMM problem
    imageResolution = size(scene_matrix);
    N_pixels = numel(scene_matrix);
    
    % Normalize true scene for metrics
    v_true_vec_norm = scene_matrix(:) / max(scene_matrix(:));
    
    % Setup operators
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
    
    % Pre-allocate variables
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    % PCG function
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    % Pre-allocate tracking arrays
    PSNR_admm_iters = zeros([config.admm_max_iter 1]);
    residuals_admm_iters = zeros([config.admm_max_iter 1]);
    
    % ADMM iterations
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
            if rel_change < 1e-6
                converged = true;
                fprintf('    ADMM converged at iteration %d\n', k_admm);
                break;
            end
        end
        
        % Print progress every 10 iterations
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
    
    % Plot convergence
    plot_admm_convergence(PSNR_admm_iters, residuals_admm_iters, config, output_folder);
end

function [Afun, Atfun_img, AtAfun_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_matrix, At_matrix, imageResolution)
    % Setup operators for ADMM
    
    % Forward operator
    Afun = @(x) A_matrix * x(:);
    
    % Adjoint operator (returns image)
    Atfun_img = @(y) reshape(At_matrix * y, imageResolution);
    
    % Normal operator (returns image)
    AtAfun_img = @(x) reshape(At_matrix * (A_matrix * x(:)), imageResolution);
    
    % TV operators
    [opDx_tv, opDtx_tv, opDtDx_tv] = difference_operators(imageResolution);
end

function [opDx_tv, opDtx_tv, opDtDx_tv] = difference_operators(imageResolution)
    % Create difference operators for TV regularization
    
    N_pixels = prod(imageResolution);
    
    % Forward difference operator (gradient)
    opDx_tv = @(x) reshape([diff(x, 1, 2), zeros(imageResolution(1), 1); ...
                            diff(x, 1, 1), zeros(1, imageResolution(2))], [N_pixels, 2]);
    
    % Adjoint difference operator (divergence)
    opDtx_tv = @(v) reshape([-v(1:N_pixels-imageResolution(1), 1); zeros(imageResolution(1), 1)] + ...
                            [v(N_pixels-imageResolution(1)+1:end, 1); zeros(N_pixels-imageResolution(1), 1)] + ...
                            [-v(1:imageResolution(1):end, 2); zeros(imageResolution(1), 1)] + ...
                            [v(imageResolution(1):imageResolution(1):end, 2); zeros(N_pixels-imageResolution(1), 1)], imageResolution);
    
    % Laplacian operator
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [x_new, z_new, u_new] = admm_iteration(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, Hfun_pcg_admm, config)
    % Single ADMM iteration
    
    % x-update
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    
    % z-update
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    z_new = max(0, 1 - kap ./ (sqrt(sum(v_z_upd.^2, 2)) + eps)) .* v_z_upd;
    
    % u-update
    u_new = u_old + opDx_tv(x_new) - z_new;
end

function [psnr, residual] = calculate_metrics(reconstructed_image, v_true_vec_norm, b_vector, Afun_admm, opDx_tv, lambda_tv_reg)
    % Calculate PSNR and residual
    
    % Normalize reconstruction
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    % PSNR
    MSE = mean((v_true_vec_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    % Residual
    data_fidelity = norm(Afun_admm(recon_norm(:)) - b_vector)^2;
    tv_penalty = lambda_tv_reg * sum(sqrt(sum(opDx_tv(recon_norm).^2, 2)));
    residual = data_fidelity + tv_penalty;
end

function plot_admm_convergence(PSNR_iters, residuals_iters, config, output_folder)
    % Plot ADMM convergence
    
    figure('Position', [100, 100, 1000, 400]);
    
    subplot(1, 2, 1);
    plot(1:length(PSNR_iters), PSNR_iters, 'b-', 'LineWidth', 2);
    xlabel('ADMM Iteration');
    ylabel('PSNR (dB)');
    title('ADMM Convergence: PSNR');
    grid on;
    
    subplot(1, 2, 2);
    plot(1:length(residuals_iters), residuals_iters, 'r-', 'LineWidth', 2);
    xlabel('ADMM Iteration');
    ylabel('Residual');
    title('ADMM Convergence: Residual');
    grid on;
    
    % Save plot
    plot_filename = fullfile(output_folder, 'admm_convergence.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  ADMM convergence plot saved: %s\n', plot_filename);
end

function plot_reconstruction_comparison(scene_matrix, reconstructed_image, X_mesh, Z_mesh, config, output_folder)
    % Plot reconstruction comparison
    
    figure('Position', [100, 100, 1200, 500]);
    
    % Convert to mm for display
    X_mm = X_mesh * 1000;
    Z_mm = Z_mesh * 1000;
    
    % Normalize images
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    % True scene
    subplot(1, 3, 1);
    imagesc(X_mm(1,:), Z_mm(:,1), scene_norm);
    colormap(gray);
    colorbar;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    title('True Scene');
    
    % Reconstructed image
    subplot(1, 3, 2);
    imagesc(X_mm(1,:), Z_mm(:,1), recon_norm);
    colormap(gray);
    colorbar;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    title('Reconstructed Image');
    
    % Difference
    subplot(1, 3, 3);
    diff_image = abs(scene_norm - recon_norm);
    imagesc(X_mm(1,:), Z_mm(:,1), diff_image);
    colormap(gray);
    colorbar;
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    title('Absolute Difference');
    
    % Calculate metrics
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    scene_vec = scene_norm(:);
    recon_vec = recon_norm(:);
    correlation = (scene_vec - mean(scene_vec))' * (recon_vec - mean(recon_vec)) / ...
        (sqrt(sum((scene_vec - mean(scene_vec)).^2)) * sqrt(sum((recon_vec - mean(recon_vec)).^2)));
    
    sgtitle(sprintf('Reconstruction Comparison (PSNR: %.2f dB, Correlation: %.4f)', psnr, correlation), 'FontSize', 14);
    
    % Save plot
    plot_filename = fullfile(output_folder, 'reconstruction_comparison.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  Reconstruction comparison plot saved: %s\n', plot_filename);
end 