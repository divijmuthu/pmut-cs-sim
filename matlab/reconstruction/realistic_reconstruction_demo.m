%% realistic_reconstruction_demo.m - Complete Realistic pMUT Reconstruction Pipeline
% Combines realistic H matrix generation (sweepScriptV29_realistic.m) 
% with full ADMM reconstruction (TheoryTest2.m)
% Uses best parameters from V29 sweep for single demonstration

clear; clc; close all;

%% ===== REALISTIC CONFIGURATION (Best from V29) =====
fprintf('=== REALISTIC pMUT RECONSTRUCTION DEMO ===\n');
fprintf('Using best parameters from V29 sweep\n\n');

% Base configuration (best from V29 sweep)
config = struct();
config.c = 343;                    % Speed of sound (m/s)
config.fs = 1e6;                   % Sampling frequency (Hz)
config.pmut_width_m = 0.020;       % pMUT width (m)
config.tx_pool_width_m = 0.200;    % Transmitter pool width (m) - BEST from V27
config.grid_width_m = 0.150;       % Imaging grid width (m)
config.target_distance_m = 0.150;  % Target distance (m)
config.grid_depth_range_m = 0.100; % Grid depth range (m)
config.grid_step_m = 0.005;        % Grid step (m) - HIGHER RESOLUTION for fine details
config.num_acquisitions = 20;      % BEST from V28_final
config.excitation_amplitude = 1e15; % Excitation amplitude

% REALISTIC pMUT PARAMETERS (from experimental data)
config.pmut_resonance_freq = 57700; % 57.7 kHz average resonance
config.pmut_bandwidth = 2520;      % 2.52 kHz standard deviation
config.impulse_duration_us = 10;   % Short impulse excitation (10 μs)

% BEST PERFORMING PARAMETERS from V29 sweep
config.num_active_tx = 5;          % Best from V28_final
config.max_delay_rand_us = 500;    % Best from V28_final
config.apodization_mode = 'uniform'; % Best from V28_final
config.frequency_offset_hz = 0;    % No frequency offset
config.num_acquisitions = 20;      % Best acquisition count

% ADMM Reconstruction Parameters (from TheoryTest2)
config.numItersADMM = 50;          % ADMM iterations
config.rho_admm = 6.73;            % Optimized ADMM penalty
config.lambda_tv_reg = 0.7438;     % Optimized TV regularization
config.admm_tol = 1.2e-5;          % Optimized tolerance
config.admm_max_iter = 50;         % Fixed max iterations
config.pcg_max_iter = 30;          % Reduced PCG iterations
config.pcg_tol = 1e-8;             % Slightly relaxed PCG tolerance
config.target_SNR_db = 35;         % Target SNR for noise

%% ===== OUTPUT SETUP =====
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('realistic_reconstruction_output', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('Output folder: %s\n\n', output_folder);

%% ===== STEP 1: REALISTIC H MATRIX GENERATION =====
fprintf('=== STEP 1: REALISTIC H MATRIX GENERATION ===\n');
fprintf('Parameters: tx=%d, delay=%dμs, apod=%s, grid=%.3fm, pool=%.3fm\n', ...
    config.num_active_tx, config.max_delay_rand_us, config.apodization_mode, ...
    config.grid_step_m, config.tx_pool_width_m);
fprintf('Impulse=%dμs, foff=%dHz, acq=%d\n\n', ...
    config.impulse_duration_us, config.frequency_offset_hz, config.num_acquisitions);

tic;
[H, coherence_matrix] = generate_h_matrix_realistic_pmut(config);
h_generation_time = toc;

fprintf('H matrix generation completed in %.2f seconds\n', h_generation_time);
fprintf('H matrix size: %d x %d\n', size(H));
fprintf('Non-zero elements: %d (%.2f%% sparsity)\n', nnz(H), (1-nnz(H)/numel(H))*100);

% Run diagnostics
fprintf('\nRunning H matrix diagnostics...\n');
metrics = run_realistic_diagnostics(H, coherence_matrix, config);

fprintf('Max Coherence: %.4f (target: <0.85)\n', metrics.max_coherence);
fprintf('Mean Coherence: %.4f\n', metrics.mean_coherence);
fprintf('Condition Number: %.2f (target: <100)\n', metrics.condition_number);
fprintf('RIP Proxy: %.4f (target: <5)\n', metrics.rip_proxy);

% Generate coherence plot
generate_realistic_coherence_plot(coherence_matrix, 'realistic_demo', output_folder, metrics);

%% ===== STEP 2: IMAGING GRID AND SCENE CREATION =====
fprintf('\n=== STEP 2: IMAGING GRID AND SCENE CREATION ===\n');

% Create imaging grid (matching H matrix dimensions)
x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
N_pixels = size(X_mesh, 1) * size(X_mesh, 2);

fprintf('Imaging grid: %d x %d = %d pixels\n', size(X_mesh, 1), size(X_mesh, 2), N_pixels);
fprintf('Grid range: X=[%.3f, %.3f]m, Z=[%.3f, %.3f]m\n', ...
    min(x_coords_img), max(x_coords_img), min(z_coords_img), max(z_coords_img));

% Create synthetic scene (similar to TheoryTest2 but adapted for X-Z grid)
scene_matrix = create_realistic_scene(X_mesh, Z_mesh, config);
v_true_vector = scene_matrix(:);

fprintf('Scene created with %d targets\n', sum(v_true_vector > 0.1));

% Plot true scene
plot_true_scene(scene_matrix, X_mesh, Z_mesh, config, output_folder);

%% ===== STEP 3: FORWARD SIMULATION =====
fprintf('\n=== STEP 3: FORWARD SIMULATION ===\n');

% Generate measurements: b = H * v_true + noise
fprintf('Generating measurements...\n');
Hv_signal = H * v_true_vector;

% Add noise
signal_power = mean(Hv_signal(:).^2);
target_SNR_linear = 10^(config.target_SNR_db / 10);
noise_variance = signal_power / target_SNR_linear;
noise_sigma = sqrt(noise_variance);
noise = noise_sigma * randn(size(Hv_signal));
b_vector = Hv_signal + noise;

% Calculate actual SNR
actual_SNR_db = 10 * log10(signal_power / mean(noise(:).^2));
fprintf('Target SNR: %.1f dB, Actual SNR: %.1f dB\n', config.target_SNR_db, actual_SNR_db);

% Plot measurement signals
plot_measurement_signals(Hv_signal, noise, b_vector, config, output_folder);

%% ===== STEP 4: ADMM RECONSTRUCTION =====
fprintf('\n=== STEP 4: ADMM RECONSTRUCTION ===\n');
fprintf('ADMM parameters: iterations=%d, rho=%.3f, lambda_TV=%.4f\n', ...
    config.admm_max_iter, config.rho_admm, config.lambda_tv_reg);

% Run ADMM reconstruction
tic;
reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, output_folder);
reconstruction_time = toc;

fprintf('Reconstruction completed in %.2f seconds\n', reconstruction_time);

%% ===== STEP 5: RESULTS VISUALIZATION =====
fprintf('\n=== STEP 5: RESULTS VISUALIZATION ===\n');

% Plot comparison
plot_reconstruction_comparison(scene_matrix, reconstructed_image, X_mesh, Z_mesh, config, output_folder);

% Calculate reconstruction metrics
calculate_reconstruction_metrics(scene_matrix, reconstructed_image, config, output_folder);

fprintf('\n=== REALISTIC RECONSTRUCTION COMPLETE ===\n');
fprintf('Results saved to: %s\n', output_folder);

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
        num_active_tx = config.num_active_tx;
        
        fprintf('    Acquisition %d/%d: Using %d transmitters...', r_acq, config.num_acquisitions, num_active_tx);
        
        % Generate REALISTIC pMUT excitation (impulse at resonant frequency)
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        
        % REALISTIC: Each pMUT has slightly different resonant frequency (from experimental data)
        individual_resonances = config.pmut_resonance_freq + config.frequency_offset_hz + ...
            (rand(1, num_active_tx) - 0.5) * config.pmut_bandwidth;
        
        % REALISTIC: Generate impulse excitation for each pMUT
        impulse_duration_samples = round(config.impulse_duration_us * fs / 1e6);
        max_len = impulse_duration_samples;
        
        % Setup apodization
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random')
            apod_weights = rand(1, num_active_tx);
        end
        
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
        
        % Setup delays
        full_delay_vector = zeros(1, vgrid_total_elements);
        if strcmp(config.apodization_mode, 'uniform')
            delays_us = zeros(1, num_active_tx);
        else
            delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        end
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

function metrics = run_realistic_diagnostics(H, coherence_matrix, config)
    % Realistic diagnostics with mutual coherence analysis
    metrics = struct();
    
    fprintf('    Computing basic metrics...\n');
    
    % Basic metrics
    metrics.sparsity = (1 - nnz(H) / numel(H)) * 100;
    
    % Column analysis
    col_norms = vecnorm(H, 2, 1);
    metrics.num_nonzero_cols = sum(col_norms > 1e-20);
    metrics.total_cols = size(H, 2);
    metrics.column_coverage = metrics.num_nonzero_cols / metrics.total_cols;
    
    fprintf('    Computing coherence...\n');
    
    % Coherence analysis
    if isempty(coherence_matrix)
        metrics.max_coherence = 0;
        metrics.mean_coherence = 0;
    else
        metrics.max_coherence = full(max(coherence_matrix(:)));
        metrics.mean_coherence = full(mean(coherence_matrix(:)));
        if isempty(metrics.max_coherence), metrics.max_coherence = 0; end
        if isempty(metrics.mean_coherence), metrics.mean_coherence = 0; end
    end
    
    fprintf('    Computing advanced metrics...\n');
    
    % Advanced metrics
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    
    if size(Hn, 2) > 1
        metrics.condition_number = cond(full(Hn));
        
        % RIP proxy
        K = min(10, size(Hn, 2));
        rip_values = zeros(3, 1);
        for i = 1:3
            idx = randperm(size(Hn, 2), K);
            s = svd(full(Hn(:, idx)));
            rip_values(i) = max(s) / min(s);
        end
        metrics.rip_proxy = mean(rip_values);
        
        % Energy concentration
        [~, S, ~] = svd(full(Hn), 'econ');
        S = diag(S);
        total_energy = sum(S.^2);
        top_energy = sum(S(1:min(10, length(S))).^2);
        metrics.energy_concentration = (top_energy / total_energy) * 100;
        
        % Rank
        metrics.rank = rank(full(Hn));
    else
        metrics.condition_number = 0;
        metrics.rip_proxy = 0;
        metrics.energy_concentration = 0;
        metrics.rank = 0;
    end
end

function generate_realistic_coherence_plot(coherence_matrix, run_name, output_folder, metrics)
    % Generate realistic coherence plot
    if isempty(coherence_matrix)
        fprintf('    Skipping coherence plot (no valid coherence matrix)\n');
        return;
    end
    
    fprintf('    Generating realistic coherence plot...\n');
    
    figure('Visible', 'off');
    
    % Create subplot layout
    subplot(2, 2, 1);
    imagesc(coherence_matrix);
    colorbar;
    title(sprintf('Realistic pMUT Coherence Matrix\nMax: %.4f, Mean: %.4f', metrics.max_coherence, metrics.mean_coherence));
    xlabel('Column Index');
    ylabel('Column Index');
    
    subplot(2, 2, 2);
    coherence_vals = coherence_matrix(:);
    coherence_vals = coherence_vals(coherence_vals > 0); % Remove zeros
    histogram(coherence_vals, 50, 'Normalization', 'probability');
    title('Realistic pMUT Coherence Distribution');
    xlabel('Coherence Value');
    ylabel('Probability');
    grid on;
    
    subplot(2, 2, 3);
    sorted_coherence = sort(coherence_vals, 'descend');
    plot(sorted_coherence(1:min(100, length(sorted_coherence))), 'b-', 'LineWidth', 2);
    title('Top 100 Coherence Values (Realistic pMUT)');
    xlabel('Rank');
    ylabel('Coherence Value');
    grid on;
    
    subplot(2, 2, 4);
    % Show coherence vs distance
    [rows, cols] = size(coherence_matrix);
    distances = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            distances(i, j) = abs(i - j);
        end
    end
    scatter(distances(:), coherence_matrix(:), 10, coherence_matrix(:), 'filled');
    colorbar;
    title('Realistic pMUT Coherence vs Column Distance');
    xlabel('Column Distance');
    ylabel('Coherence Value');
    grid on;
    
    % Save plot
    plot_filename = fullfile(output_folder, sprintf('realistic_coherence_plot_%s.png', run_name));
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('    Realistic coherence plot saved: %s\n', plot_filename);
end

function scene_matrix = create_realistic_scene(X_mesh, Z_mesh, config)
    % Create complex synthetic scene with 3x3 grid of white dots on black background
    
    scene_matrix = zeros(size(X_mesh));
    
    % Convert to mm for target placement
    X_mm = X_mesh * 1000;
    Z_mm = Z_mesh * 1000;
    
    % Debug: Print grid coordinate ranges
    fprintf('    X grid range: [%.1f, %.1f] mm\n', min(X_mm(:)), max(X_mm(:)));
    fprintf('    Z grid range: [%.1f, %.1f] mm\n', min(Z_mm(:)), max(Z_mm(:)));
    fprintf('    X grid step: %.1f mm\n', (X_mm(1,2) - X_mm(1,1)));
    fprintf('    Z grid step: %.1f mm\n', (Z_mm(2,1) - Z_mm(1,1)));
    
    % Define 3x3 grid of target positions (9 dots total)
    % Grid spacing: 20mm between dots (closer spacing for more challenging reconstruction)
    grid_spacing_mm = 20;
    grid_offset_x_mm = 0;   % Center the grid
    grid_offset_z_mm = 150; % Position in middle of Z range (100-200mm)
    
    target_positions = [];
    for row = 1:3
        for col = 1:3
            x_pos_mm = grid_offset_x_mm + (col - 2) * grid_spacing_mm;  % Center at 0
            z_pos_mm = grid_offset_z_mm + (row - 2) * grid_spacing_mm;  % Center at 50mm
            target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
        end
    end
    
    % Target size (smaller dots for fine details and potential mixing)
    target_size_mm = 4;  % 4mm dots for fine details and challenging reconstruction
    grid_step_mm = config.grid_step_m * 1000;
    target_radius_pixels = round(target_size_mm / (2 * grid_step_mm));
    
    % Ensure minimum size
    if target_radius_pixels < 1
        target_radius_pixels = 1;
    end
    
    for i = 1:size(target_positions, 1)
        x_pos_mm = target_positions(i, 1);
        z_pos_mm = target_positions(i, 2);
        amplitude = target_positions(i, 3);
        
        % Find closest grid points
        [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
        [~, iz_scene] = min(abs(Z_mm(:,1) - z_pos_mm));
        
        % Debug: Print the actual grid coordinates
        fprintf('    Target %d: desired (%.1f, %.1f) mm -> grid (%.1f, %.1f) mm -> indices (%d, %d)\n', ...
            i, x_pos_mm, z_pos_mm, X_mm(1,ix_scene), Z_mm(iz_scene,1), ix_scene, iz_scene);
        
            % Create target dot
    for dz = -target_radius_pixels:target_radius_pixels
        for dx = -target_radius_pixels:target_radius_pixels
            ix_target = ix_scene + dx;
            iz_target = iz_scene + dz;
            
            % Check bounds
            if ix_target > 0 && ix_target <= size(X_mesh, 2) && ...
               iz_target > 0 && iz_target <= size(X_mesh, 1)
                scene_matrix(iz_target, ix_target) = amplitude;
            end
        end
    end
    
    % Debug: Print actual grid indices for this target
    fprintf('    Target %d placed at grid indices: (%d, %d) for position (%.1f, %.1f) mm\n', ...
        i, ix_scene, iz_scene, x_pos_mm, z_pos_mm);
    end
    
    fprintf('  Complex scene created with 3x3 grid of 9 targets:\n');
    fprintf('    Grid spacing: %.1f mm\n', grid_spacing_mm);
    fprintf('    Target size: %.1f mm\n', target_size_mm);
    fprintf('    Grid center: (%.1f, %.1f) mm\n', grid_offset_x_mm, grid_offset_z_mm);
    fprintf('    Total targets: %d\n', size(target_positions, 1));
    
    % Print target positions
    for i = 1:size(target_positions, 1)
        fprintf('    Target %d: (%.1f, %.1f) mm, amplitude %.1f\n', ...
            i, target_positions(i, 1), target_positions(i, 2), target_positions(i, 3));
    end
end

function plot_true_scene(scene_matrix, X_mesh, Z_mesh, config, output_folder)
    % Plot the true scene
    figure('Visible', 'off');
    
    % Debug: Print scene matrix info
    fprintf('    Scene matrix size: %d x %d\n', size(scene_matrix));
    fprintf('    Scene matrix range: [%.3f, %.3f]\n', min(scene_matrix(:)), max(scene_matrix(:)));
    fprintf('    Number of non-zero pixels: %d\n', nnz(scene_matrix));
    
    % Main scene plot with proper colormap
    subplot(1, 2, 1);
    imagesc(X_mesh(1,:)*1000, Z_mesh(:,1)*1000, scene_matrix);
    axis image;
    colormap gray;
    colorbar;
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    title('True Scene (Realistic pMUT Demo)');
    set(gca, 'YDir', 'normal');
    
    % Vector plot to show individual targets
    subplot(1, 2, 2);
    stem(scene_matrix(:));
    axis square tight;
    title('True Scene (vector)');
    xlabel('Pixel Index');
    ylabel('Amplitude');
    grid on;
    
    % Save plot
    plot_filename = fullfile(output_folder, 'true_scene.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  True scene plot saved: %s\n', plot_filename);
end

function plot_measurement_signals(Hv_signal, noise, b_vector, config, output_folder)
    % Plot measurement signals
    figure('Visible', 'off');
    
    subplot(3, 1, 1);
    plot(Hv_signal(1:min(1000, length(Hv_signal))));
    title('Clean Signal (H * v_{true})');
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on;
    
    subplot(3, 1, 2);
    plot(noise(1:min(1000, length(noise))));
    title('Noise');
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on;
    
    subplot(3, 1, 3);
    plot(b_vector(1:min(1000, length(b_vector))));
    title('Measured Signal (H * v_{true} + noise)');
    xlabel('Sample Index');
    ylabel('Amplitude');
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
    
    % Pre-compute matrices
    A_matrix = H;
    I_true_matrix = scene_matrix;
    v_true_vec_norm = scene_matrix(:);
    if max(abs(v_true_vec_norm)) > 0
        v_true_vec_norm = v_true_vec_norm ./ max(abs(v_true_vec_norm));
    end
    imageResolution = size(I_true_matrix);
    
    % Check dimensions
    if size(A_matrix, 2) ~= length(v_true_vec_norm)
        fprintf('  ERROR: H matrix columns (%d) != scene vector length (%d)\n', size(A_matrix, 2), length(v_true_vec_norm));
        reconstructed_image = zeros(imageResolution);
        return;
    end
    
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
    
    % Setup visualization
    figure('Visible', 'off');
    set(gcf, 'Position', [200, 200, 1000, 700], 'Color', 'w');
    
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
    
    % Plot convergence
    plot_convergence(PSNR_admm_iters, residuals_admm_iters, config, output_folder);
    
    close(gcf);
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    % Setup operators for ADMM
    
    % Matrix operators
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(b) reshape(At_admm * b, imageResolution);
    AtAfun_admm_img = @(x) reshape(At_admm * (A_admm * x(:)), imageResolution);
    
    % TV operators
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
    % Single ADMM iteration
    
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

function plot_convergence(PSNR_admm_iters, residuals_admm_iters, config, output_folder)
    % Plot convergence curves
    
    figure('Visible', 'off');
    
    subplot(2, 1, 1);
    plot(1:length(PSNR_admm_iters), PSNR_admm_iters, 'b-', 'LineWidth', 2);
    title('ADMM Convergence: PSNR vs Iteration');
    xlabel('Iteration');
    ylabel('PSNR (dB)');
    grid on;
    
    subplot(2, 1, 2);
    plot(1:length(residuals_admm_iters), residuals_admm_iters, 'r-', 'LineWidth', 2);
    title('ADMM Convergence: Residual vs Iteration');
    xlabel('Iteration');
    ylabel('Residual');
    grid on;
    
    % Save plot
    plot_filename = fullfile(output_folder, 'admm_convergence.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  ADMM convergence plot saved: %s\n', plot_filename);
end

function plot_reconstruction_comparison(scene_matrix, reconstructed_image, X_mesh, Z_mesh, config, output_folder)
    % Plot comparison between true and reconstructed images
    
    figure('Visible', 'off');
    
    subplot(2, 2, 1);
    imagesc(X_mesh(1,:)*1000, Z_mesh(:,1)*1000, scene_matrix);
    colorbar;
    title('True Scene');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    axis equal tight;
    
    subplot(2, 2, 2);
    imagesc(X_mesh(1,:)*1000, Z_mesh(:,1)*1000, reconstructed_image);
    colorbar;
    title('Reconstructed Image');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    axis equal tight;
    
    subplot(2, 2, 3);
    % Difference image
    diff_image = abs(reconstructed_image - scene_matrix);
    imagesc(X_mesh(1,:)*1000, Z_mesh(:,1)*1000, diff_image);
    colorbar;
    title('Absolute Difference');
    xlabel('X Position (mm)');
    ylabel('Z Position (mm)');
    axis equal tight;
    
    subplot(2, 2, 4);
    % Cross-section comparison
    [~, max_idx] = max(scene_matrix(:));
    [iz, ix] = ind2sub(size(scene_matrix), max_idx);
    
    plot(X_mesh(1,:)*1000, scene_matrix(iz,:), 'b-', 'LineWidth', 2, 'DisplayName', 'True');
    hold on;
    plot(X_mesh(1,:)*1000, reconstructed_image(iz,:), 'r--', 'LineWidth', 2, 'DisplayName', 'Reconstructed');
    title(sprintf('Cross-section at Z=%.1f mm', Z_mesh(iz,1)*1000));
    xlabel('X Position (mm)');
    ylabel('Amplitude');
    legend('Location', 'best');
    grid on;
    
    % Save plot
    plot_filename = fullfile(output_folder, 'reconstruction_comparison.png');
    saveas(gcf, plot_filename);
    close(gcf);
    
    fprintf('  Reconstruction comparison plot saved: %s\n', plot_filename);
end

function calculate_reconstruction_metrics(scene_matrix, reconstructed_image, config, output_folder)
    % Calculate and display reconstruction metrics
    
    % Normalize images
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    % PSNR
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    PSNR = 10 * log10(1 / MSE);
    
    % SSIM (if available)
    try
        SSIM_val = ssim(recon_norm, scene_norm);
    catch
        SSIM_val = NaN;
    end
    
    % Correlation (manual calculation to avoid toolbox dependency)
    scene_vec = scene_norm(:);
    recon_vec = recon_norm(:);
    correlation = (scene_vec - mean(scene_vec))' * (recon_vec - mean(recon_vec)) / ...
        (sqrt(sum((scene_vec - mean(scene_vec)).^2)) * sqrt(sum((recon_vec - mean(recon_vec)).^2)));
    
    % Print metrics
    fprintf('\n=== RECONSTRUCTION METRICS ===\n');
    fprintf('PSNR: %.2f dB\n', PSNR);
    if ~isnan(SSIM_val)
        fprintf('SSIM: %.4f\n', SSIM_val);
    end
    fprintf('Correlation: %.4f\n', correlation);
    
    % Save metrics to file
    metrics_file = fullfile(output_folder, 'reconstruction_metrics.txt');
    fid = fopen(metrics_file, 'w');
    fprintf(fid, 'Reconstruction Metrics\n');
    fprintf(fid, '=====================\n\n');
    fprintf(fid, 'PSNR: %.2f dB\n', PSNR);
    if ~isnan(SSIM_val)
        fprintf(fid, 'SSIM: %.4f\n', SSIM_val);
    end
    fprintf(fid, 'Correlation: %.4f\n', correlation);
    fprintf(fid, '\nConfiguration:\n');
    fprintf(fid, 'Number of transmitters: %d\n', config.num_active_tx);
    fprintf(fid, 'Delay range: %d μs\n', config.max_delay_rand_us);
    fprintf(fid, 'Apodization: %s\n', config.apodization_mode);
    fprintf(fid, 'Grid step: %.3f m\n', config.grid_step_m);
    fprintf(fid, 'Impulse duration: %d μs\n', config.impulse_duration_us);
    fprintf(fid, 'Frequency offset: %d Hz\n', config.frequency_offset_hz);
    fprintf(fid, 'Number of acquisitions: %d\n', config.num_acquisitions);
    fclose(fid);
    
    fprintf('  Reconstruction metrics saved to: %s\n', metrics_file);
end 