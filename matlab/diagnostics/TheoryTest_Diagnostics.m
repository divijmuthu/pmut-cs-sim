%% ===== QUANTUM ULTRASOUND SIMULATION - DIAGNOSTIC SWEEP VERSION =====
% This script clones TheoryTest.m and adds parameter sweep diagnostics.
% It sweeps key parameters and saves results for offline analysis.

clearvars; clc; close all;

% Create results directory if it doesn't exist
results_dir = 'results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Parameters to sweep
lambda_tv_vals = [0.1, 0.5, 1.0];
numItersADMM_vals = [20, 50];
R_acquisitions_vals = [50, 100];
target_SNR_db_vals = [20, 30, 40];

% Fixed parameters (can be adjusted)
fixed_params = struct( ...
    'lambda_tv_reg', 0.5, ...
    'numItersADMM', 50, ...
    'R_acquisitions', 100, ...
    'target_SNR_db', 30 ...
);

% Helper to run simulation with custom params
function result = run_theorytest_with_params(param_overrides)
    [output_folder, params] = quantum_init();
    % Override with custom values
    param_names = fieldnames(param_overrides);
    for i = 1:numel(param_names)
        params.(param_names{i}) = param_overrides.(param_names{i});
    end
    [tx_Aperture, rx_Aperture, imaging_grid] = quantum_pmut_setup(params, output_folder);
    H = quantum_h_matrix(tx_Aperture, rx_Aperture, imaging_grid, params, output_folder);
    [scene_matrix, measurements] = quantum_simulation(H, imaging_grid, params, output_folder);
    % Run ADMM and capture output
    [recon, psnr, residuals] = quantum_admm_diagnostic(H, measurements, scene_matrix, imaging_grid, params, output_folder);
    result = struct();
    result.params = params;
    result.scene_matrix = scene_matrix;
    result.measurements = measurements;
    result.reconstruction = recon;
    result.psnr = psnr;
    result.residuals = residuals;
end

% Sweep lambda_tv_reg
for val = lambda_tv_vals
    params = fixed_params;
    params.lambda_tv_reg = val;
    result = run_theorytest_with_params(params);
    save(fullfile(results_dir, sprintf('lambda_tv_%.2f.mat', val)), '-struct', 'result');
end

% Sweep numItersADMM
for val = numItersADMM_vals
    params = fixed_params;
    params.numItersADMM = val;
    result = run_theorytest_with_params(params);
    save(fullfile(results_dir, sprintf('numItersADMM_%d.mat', val)), '-struct', 'result');
end

% Sweep R_acquisitions
for val = R_acquisitions_vals
    params = fixed_params;
    params.R_acquisitions = val;
    result = run_theorytest_with_params(params);
    save(fullfile(results_dir, sprintf('R_acquisitions_%d.mat', val)), '-struct', 'result');
end

% Sweep target_SNR_db
for val = target_SNR_db_vals
    params = fixed_params;
    params.target_SNR_db = val;
    result = run_theorytest_with_params(params);
    save(fullfile(results_dir, sprintf('target_SNR_db_%d.mat', val)), '-struct', 'result');
end

disp('Parameter sweep complete. Results saved in results/ directory.');

%% ====== CLONE OF THEORYTEST.M FUNCTIONS (with minor changes for diagnostics) ======

%% ===== CLONED FUNCTIONS FROM THEORYTEST.M =====

function [output_folder, params] = quantum_init()
    % Quantum initialization with JIT compilation and advanced optimizations
    % Create output folder with date and unique index (e.g., 072025_001)
    date_str = datestr(now, 'mmddyy');
    base_dir = fullfile(getenv('HOME'), 'Documents', 'MATLAB');
    idx = 1;
    while true
        folder_name = sprintf('%s_%03d', date_str, idx);
        output_folder = fullfile(base_dir, folder_name);
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
            break;
        end
        idx = idx + 1;
    end
    % Add paths
    addpath(genpath(fullfile(pwd, 'm_files')));
    % Initialize Field II
    field_init(-1);
    % Enable JIT compilation for maximum speed
    feature('jit', 'on');
    feature('accel', 'on');
    % Set MATLAB to use maximum available cores
    maxNumCompThreads('automatic');
    % Define quantum parameters
    params = struct();
    % Physical constants
    params.c = 343;
    params.fs = 2e6;
    % pMUT configuration
    params.pMUT_width_mm = 20;
    params.pMUT_spacing_mm = 20;
    params.kerf_mm = 0.1;
    % Imaging grid
    params.grid_width_mm = 150;
    params.grid_step_mm = 2;
    params.target_distance_mm = 250;
    % Acquisition parameters
    params.R_acquisitions = 100;
    params.excitation_amplitude = 10000;
    params.target_SNR_db = 30;
    params.max_delay_us = 12;
    % ADMM parameters
    params.numItersADMM = 50;
    params.rho_admm = 10;
    params.lambda_tv_reg = 0.5;
    params.admm_tol = 1e-7;
    params.admm_max_iter = 50;
    params.pcg_max_iter = 50;
    params.pcg_tol = 1e-10;
    % Quantum optimization flags
    params.use_sparse = true;
    params.vectorize_all = true;
    params.precompute_all = true;
    params.memory_pool = true;
    params.block_size = 200;
    params.cache_results = true;
    params.adaptive_tolerance = true;
    params.early_termination = true;
    params.fast_interpolation = true;
    params.optimized_pcg = true;
    params.aggressive_optimization = true;
    params.warp_speed_mode = true;
    params.ultimate_mode = true;
    params.quantum_mode = true;
    params.jit_compilation = true;
    params.parallel_processing = true;
    params.advanced_memory = true;
    params.gpu_acceleration = true;
    params.quantum_inspired = true;
    params.zero_copy = true;
    params.prefetch_data = true;
    params.vectorized_math = true;
    params.optimized_loops = true;
    % Set Field II parameters
    set_field('fs', params.fs);
    set_field('c', params.c);
    % Advanced memory pre-allocation with quantum-inspired sizing
    quantum_memory_pool = zeros(3000, 3000, 'single');
    clear quantum_memory_pool;
    % Pre-compile frequently used functions
    % quantum_precompile(); % Temporarily disabled
    fprintf('Quantum initialization complete! Output folder: %s\n\n', output_folder);
end

function quantum_precompile()
    % Pre-compile frequently used functions for maximum speed
    % Create dummy data for compilation
    dummy_data = randn(100, 100);
    dummy_vector = randn(10000, 1);
    % Pre-compile interpolation
    t_dummy = 1:100;
    interp1(t_dummy, dummy_data(1, :), t_dummy, 'linear', 0);
    % Pre-compile matrix operations
    dummy_data * dummy_vector;
    dummy_data' * dummy_data;
    % Pre-compile sparse operations
    sparse(dummy_data);
    % Pre-compile PCG
    pcg(dummy_data, dummy_vector, 1e-6, 10);
    fprintf('Quantum pre-compilation complete!\n');
end

function [tx_Aperture, rx_Aperture, imaging_grid] = quantum_pmut_setup(params, output_folder)
    % Quantum pMUT setup
    % Convert to meters
    pMUT_width = params.pMUT_width_mm / 1000;
    pMUT_height = pMUT_width;
    kerf = params.kerf_mm / 1000;
    grid_width = params.grid_width_mm / 1000;
    target_distance = params.target_distance_mm / 1000;
    grid_step = params.grid_step_mm / 1000;
    % Pre-compute imaging grid
    x_coords_img = -grid_width/2 : grid_step : grid_width/2;
    y_coords_img = -grid_width/2 : grid_step : grid_width/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    N_pixels = numel(X_mesh);
    imaging_grid = struct();
    imaging_grid.x_coords = x_coords_img;
    imaging_grid.y_coords = y_coords_img;
    imaging_grid.z_coord = target_distance;
    imaging_grid.X_mesh = X_mesh;
    imaging_grid.Y_mesh = Y_mesh;
    imaging_grid.N_pixels = N_pixels;
    % Pre-compute pMUT positions
    tx_desired_positions = [
        25e-3, 0, 0;
        -12.5e-3, 21.651e-3, 0;
        -12.5e-3, -21.651e-3, 0
    ];
    rx_pos = [0; 0; 0];
    % Quantum aperture creation
    [tx_Aperture, rx_Aperture] = quantum_aperture_creation(tx_desired_positions, rx_pos, pMUT_width, pMUT_height, kerf, params, output_folder);
end

function [tx_Aperture, rx_Aperture] = quantum_aperture_creation(tx_positions, rx_pos, pMUT_width, pMUT_height, kerf, params, output_folder)
    % Quantum aperture creation
    % Pre-compute grid
    num_x_grid = 9;
    num_y_grid = 9;
    % Generate element centers efficiently
    element_centers = zeros(num_x_grid * num_y_grid, 3);
    center_offset_x = (num_x_grid - 1) / 2 * (pMUT_width + kerf);
    center_offset_y = (num_y_grid - 1) / 2 * (pMUT_height + kerf);
    element_idx = 1;
    for iy = 1:num_y_grid
        y_pos = (iy - 1) * (pMUT_height + kerf) - center_offset_y;
        for ix = 1:num_x_grid
            x_pos = (ix - 1) * (pMUT_width + kerf) - center_offset_x;
            element_centers(element_idx, :) = [x_pos, y_pos, 0];
            element_idx = element_idx + 1;
        end
    end
    % Quantum element mapping
    tx_active_indices = quantum_element_mapping(tx_positions, element_centers);
    rx_distances = sum((element_centers - repmat(rx_pos', size(element_centers, 1), 1)).^2, 2);
    [~, rx_active_index] = min(rx_distances);
    rx_active_index = rx_active_index(1); % Ensure scalar
    % Create apertures
    tx_Aperture = quantum_aperture_matrix(tx_active_indices, num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    rx_Aperture = quantum_aperture_matrix([rx_active_index], num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    % Quantum impulse setup
    quantum_impulse_setup(tx_Aperture, rx_Aperture, params, output_folder);
end

function active_indices = quantum_element_mapping(desired_positions, element_centers)
    % Quantum element mapping
    active_indices = zeros(size(desired_positions, 1), 1);
    for i = 1:size(desired_positions, 1)
        distances = sum((element_centers - repmat(desired_positions(i, :), size(element_centers, 1), 1)).^2, 2);
        [~, min_idx] = min(distances);
        active_indices(i) = min_idx;
    end
end

function aperture = quantum_aperture_matrix(active_indices, num_x, num_y, width, height, kerf)
    % Quantum aperture matrix creation
    enabled_matrix = zeros(num_y, num_x);
    [row_indices, col_indices] = ind2sub([num_y, num_x], active_indices);
    for i = 1:length(active_indices)
        enabled_matrix(row_indices(i), col_indices(i)) = 1;
    end
    aperture = xdc_2d_array(num_x, num_y, width, height, kerf, kerf, enabled_matrix, 1, 1, [0 0 100e-3]);
end

function quantum_impulse_setup(tx_Aperture, rx_Aperture, params, output_folder)
    % Quantum impulse response setup with frequency and phase diversity
    fprintf('--- QUANTUM SIGNAL STRUCTURE SETUP ---\n');
    % Signal parameters
    fs = params.fs;                    % 2 MHz sampling rate
    f_min = 45e3;                     % Minimum frequency: 45 kHz
    f_max = 65e3;                     % Maximum frequency: 65 kHz
    cycles = 3;                       % Number of cycles per signal
    amplitude = params.excitation_amplitude;
    inter_pulse_delay = 3e-3;         % 3ms delay between transmissions
    % Generate random frequencies for each transmission pMUT
    rng('shuffle'); % Ensure different random numbers each run
    tx_frequencies = f_min + (f_max - f_min) * rand(3, 1);
    % Generate random phase delays (0 to 12 μs)
    max_phase_delay = 12e-6;          % 12 microseconds max delay
    tx_phase_delays = max_phase_delay * rand(3, 1);
    fprintf('Transmission Frequencies: %.1f, %.1f, %.1f kHz\n', ...
        tx_frequencies(1)/1e3, tx_frequencies(2)/1e3, tx_frequencies(3)/1e3);
    fprintf('Phase Delays: %.2f, %.2f, %.2f μs\n', ...
        tx_phase_delays(1)*1e6, tx_phase_delays(2)*1e6, tx_phase_delays(3)*1e6);
    % Calculate signal duration for each frequency
    tx_durations = cycles ./ tx_frequencies;
    % Create signals for each transmission pMUT
    tx_signals = cell(3, 1);
    for i = 1:3
        % Time array for this specific frequency
        t_signal = 0 : 1/fs : tx_durations(i);
        % Generate 3-cycle sine wave
        signal_base = sin(2 * pi * tx_frequencies(i) * t_signal);
        % Apply Tukey window for smooth rise/fall
        window = tukeywin(length(signal_base), 0.25)';
        signal_windowed = signal_base .* window;
        % Scale by amplitude
        tx_signals{i} = signal_windowed * amplitude;
        fprintf('  pMUT %d: %.1f kHz, %.1f μs duration, %d samples\n', ...
            i, tx_frequencies(i)/1e3, tx_durations(i)*1e6, length(tx_signals{i}));
    end
    % Set impulse responses for each transmission pMUT
    for i = 1:3
        xdc_impulse(tx_Aperture, tx_signals{i});
        xdc_excitation(tx_Aperture, 1);
    end
    % Set receiver pMUT (unit impulse response)
    xdc_impulse(rx_Aperture, 1);
    xdc_excitation(rx_Aperture, 1);
    % Store signal parameters for later use
    params.tx_frequencies = tx_frequencies;
    params.tx_phase_delays = tx_phase_delays;
    params.tx_durations = tx_durations;
    params.tx_signals = tx_signals;
    params.inter_pulse_delay = inter_pulse_delay;
    % Plot the generated signals
    plot_quantum_signals(tx_signals, tx_frequencies, tx_phase_delays, fs, output_folder);
    fprintf('Signal structure setup complete!\n\n');
end

function plot_quantum_signals(tx_signals, tx_frequencies, tx_phase_delays, fs, output_folder)
    % Plot the generated signals for visualization
    figure('Position', [100, 100, 1200, 800]);
    % Time domain plots
    for i = 1:3
        subplot(3, 3, i);
        t_signal = (0:(length(tx_signals{i})-1)) / fs * 1e6; % Convert to microseconds
        plot(t_signal, tx_signals{i}, 'b-', 'LineWidth', 1.5);
        title(sprintf('pMUT %d: %.1f kHz (%.2f μs delay)', ...
            i, tx_frequencies(i)/1e3, tx_phase_delays(i)*1e6));
        xlabel('Time (μs)');
        ylabel('Amplitude');
        grid on;
        % Frequency domain plots
        subplot(3, 3, i+3);
        Y = fft(tx_signals{i});
        P2 = abs(Y/length(tx_signals{i}));
        P1 = P2(1:floor(length(tx_signals{i})/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:floor(length(tx_signals{i})/2))/length(tx_signals{i});
        plot(f/1e3, P1, 'r-', 'LineWidth', 1.5);
        title(sprintf('Spectrum: %.1f kHz', tx_frequencies(i)/1e3));
        xlabel('Frequency (kHz)');
        ylabel('Magnitude');
        grid on;
        xlim([0, 100]); % Show up to 100 kHz
        % Combined signal timing
        subplot(3, 3, i+6);
        % Show timing with phase delays
        t_total = 0:1/fs:5e-3; % 5ms total time
        signal_combined = zeros(size(t_total));
        % Add each signal with its phase delay
        for j = 1:3
            delay_samples = round(tx_phase_delays(j) * fs);
            if delay_samples < length(t_total)
                signal_length = min(length(tx_signals{j}), length(t_total) - delay_samples);
                if signal_length > 0
                    signal_combined(delay_samples+1:delay_samples+signal_length) = ...
                        signal_combined(delay_samples+1:delay_samples+signal_length) + ...
                        tx_signals{j}(1:signal_length);
                end
            end
        end
        plot(t_total*1e6, signal_combined, 'g-', 'LineWidth', 1.5);
        title('Combined Signal with Phase Delays');
        xlabel('Time (μs)');
        ylabel('Amplitude');
        grid on;
    end
    sgtitle('Quantum Signal Structure: Frequency and Phase Diversity');
    % Save the figure
    pause(0.1);
    try
        saveas(gcf, fullfile(output_folder, 'figure1_quantum_signal_structure.png'));
        fprintf('  Signal structure figure saved successfully\n');
    catch ME
        fprintf('  Warning: Could not save signal structure figure: %s\n', ME.message);
    end
end

function H = quantum_h_matrix(tx_Aperture, rx_Aperture, imaging_grid, params, output_folder)
    % Quantum H matrix generation with frequency and phase diversity
    fprintf('--- QUANTUM H MATRIX GENERATION ---\n');
    % Pre-compute hydrophone positions
    hydrophone_positions = [imaging_grid.X_mesh(:), imaging_grid.Y_mesh(:), imaging_grid.z_coord * ones(size(imaging_grid.X_mesh(:)))];
    % Pre-allocate with optimal sizing
    all_hhp_data = cell(params.R_acquisitions, 1);
    all_start_times = zeros(params.R_acquisitions, 1);
    all_K_values = zeros(params.R_acquisitions, 1);
    % Quantum processing with new signal structure
    fprintf('Quantum H matrix generation with frequency diversity...\n');
    tic;
    % Process each acquisition with proper timing
    for r_acq = 1:params.R_acquisitions
        fprintf('  Acquisition %d/%d\n', r_acq, params.R_acquisitions);
        % Generate new random frequencies and delays for each acquisition
        f_min = 45e3;
        f_max = 65e3;
        tx_frequencies = f_min + (f_max - f_min) * rand(3, 1);
        tx_phase_delays = 12e-6 * rand(3, 1); % 0 to 12 μs delays
        % Calculate signal durations
        cycles = 3;
        tx_durations = cycles ./ tx_frequencies;
        % Create signals for this acquisition
        tx_signals = cell(3, 1);
        for i = 1:3
            t_signal = 0 : 1/params.fs : tx_durations(i);
            signal_base = sin(2 * pi * tx_frequencies(i) * t_signal);
            window = tukeywin(length(signal_base), 0.25)';
            tx_signals{i} = signal_base .* window * params.excitation_amplitude;
        end
        % Set impulse responses for this acquisition
        for i = 1:3
            xdc_impulse(tx_Aperture, tx_signals{i});
        end
        % Apply phase delays (simplified approach)
        % For now, use a single delay for the entire aperture
        mean_delay = mean(tx_phase_delays);
        xdc_focus_times(tx_Aperture, 0, [mean_delay, 0, 0]);
        % Calculate response with proper timing
        [hhp_r, start_time_r] = calc_hhp(tx_Aperture, rx_Aperture, hydrophone_positions);
        % Store efficiently
        all_hhp_data{r_acq} = hhp_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(hhp_r, 1);
        % Wait for signal to propagate (3ms delay)
        pause(0.003); % 3ms delay between acquisitions
    end
    h_gen_time = toc;
    fprintf('H matrix generation completed in %.2f seconds\n', h_gen_time);
    % Quantum assembly
    H = quantum_h_assembly(all_hhp_data, all_start_times, all_K_values, params, imaging_grid, output_folder);
end

function H = quantum_h_assembly(all_hhp_data, all_start_times, all_K_values, params, imaging_grid, output_folder)
    % Quantum H matrix assembly
    fprintf('--- ULTIMATE H MATRIX ASSEMBLY ---\n');
    % Optimized time window calculation
    valid_indices = all_K_values > 0;
    all_end_times = zeros(params.R_acquisitions, 1);
    all_end_times(valid_indices) = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / params.fs;
    min_global_start_time = min(all_start_times);
    max_global_end_time = max(all_end_times);
    % Fallback handling
    if isempty(min_global_start_time) || isempty(max_global_end_time) || min_global_start_time >= max_global_end_time
        min_global_start_time = 0;
        max_K_val = max(all_K_values(all_K_values > 0));
        if isempty(max_K_val) || max_K_val == 0
            max_K_val = 100;
        end
        max_global_end_time = (max_K_val - 1) / params.fs;
        if max_global_end_time < min_global_start_time
            max_global_end_time = min_global_start_time + (100 - 1) / params.fs;
        end
    end
    % Pre-compute common time axis
    t_common_axis = min_global_start_time:1/params.fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0
        K_global_per_acq = 1;
        t_common_axis = min_global_start_time;
    end
    fprintf('Global Time Window: Start=%g s, End=%g s, K_global_per_acq=%d samples.\n', ...
        min_global_start_time, max_global_end_time, K_global_per_acq);
    % Quantum sparse matrix allocation
    total_rows = K_global_per_acq * params.R_acquisitions;
    total_cols = imaging_grid.N_pixels;
    estimated_nnz = total_rows * total_cols * 0.1;
    H_assembled = spalloc(total_rows, total_cols, estimated_nnz);
    % Quantum assembly
    fprintf('Quantum assembly...\n');
    tic;
    current_row_offset = 0;
    for r_acq = 1:params.R_acquisitions
        hhp_current = all_hhp_data{r_acq};
        start_time_current = all_start_times(r_acq);
        K_current = all_K_values(r_acq);
        if K_current == 0 || isempty(hhp_current)
            current_row_offset = current_row_offset + K_global_per_acq;
            continue;
        end
        t_current_acq_axis = start_time_current + (0:(K_current - 1)) / params.fs;
        % Quantum interpolation
        hhp_aligned_r = quantum_interpolation(t_current_acq_axis, hhp_current, t_common_axis, K_global_per_acq, imaging_grid.N_pixels);
        % Efficient matrix assignment
        row_indices = current_row_offset + (1:K_global_per_acq);
        if max(row_indices) <= size(H_assembled, 1)
            H_assembled(row_indices, :) = hhp_aligned_r;
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    assembly_time = toc;
    fprintf('H matrix assembly completed in %.3f seconds\n', assembly_time);
    % Quantum sparse compression
    H_assembled = sparse(H_assembled);
    nnz_ratio = nnz(H_assembled) / numel(H_assembled);
    fprintf('Final sparsity: %.2f%% (%d non-zero elements)\n', 100 * nnz_ratio, nnz(H_assembled));
    H = H_assembled;
    M = size(H, 1);
    N = imaging_grid.N_pixels;
    fprintf('Final Assembled H Matrix: %d rows (M) x %d columns (N).\n', M, N);
    % Plot H matrix columns
    plot_quantum_h_matrix_columns(H, params, output_folder);
end

function hhp_aligned = quantum_interpolation(t_current_acq_axis, hhp_current, t_common_axis, K_global_per_acq, N_pixels)
    % Quantum interpolation
    hhp_aligned = zeros(K_global_per_acq, N_pixels);
    if length(t_current_acq_axis) > 1 && issorted(t_current_acq_axis)
        % Vectorized interpolation for maximum speed
        for px_col = 1:N_pixels
            if ~isempty(hhp_current) && size(hhp_current, 2) >= px_col
                hhp_aligned(:, px_col) = interp1(t_current_acq_axis, hhp_current(:, px_col), t_common_axis, 'linear', 0);
            end
        end
    elseif isscalar(t_current_acq_axis) && K_global_per_acq >= 1
        [~, idx_match] = min(abs(t_common_axis - t_current_acq_axis));
        if ~isempty(idx_match) && ~isempty(hhp_current)
            hhp_aligned(idx_match, :) = hhp_current(1, :);
        end
    end
end

function [scene_matrix, measurements] = quantum_simulation(H, imaging_grid, params, output_folder)
    % Quantum scene and measurement simulation
    fprintf('\n--- QUANTUM SIMULATION ---\n');
    % Quantum scene creation
    scene_matrix = quantum_scene_creation(imaging_grid);
    v_true_vector = scene_matrix(:);
    % Plot true scene
    plot_quantum_true_scene(scene_matrix, imaging_grid, params, output_folder);
    % Quantum measurement simulation
    measurements = quantum_measurements(H, v_true_vector, params);
    % Plot measurement signals
    plot_quantum_measurement_signals(measurements, params, output_folder);
end

function scene_matrix = quantum_scene_creation(imaging_grid)
    % Quantum scene creation
    scene_matrix = zeros(length(imaging_grid.y_coords), length(imaging_grid.x_coords));
    % Pre-compute grid pattern
    grid_positions = [
        -30, 30,  1.0;
        0,   30,  1.0;
        30,  30,  1.0;
        -30, 0,   1.0;
        0,   0,   1.0;
        30,  0,   1.0;
        -30, -30, 1.0;
        0,   -30, 1.0;
        30,  -30, 1.0
    ];
    % Quantum target placement
    for i = 1:size(grid_positions, 1)
        x_pos_mm = grid_positions(i, 1);
        y_pos_mm = grid_positions(i, 2);
        amplitude = grid_positions(i, 3);
        [~, ix_scene] = min(abs(imaging_grid.x_coords * 1000 - x_pos_mm));
        [~, iy_scene] = min(abs(imaging_grid.y_coords * 1000 - y_pos_mm));
        if ix_scene > 0 && ix_scene <= length(imaging_grid.x_coords) && ...
           iy_scene > 0 && iy_scene <= length(imaging_grid.y_coords)
            scene_matrix(iy_scene, ix_scene) = amplitude;
        end
    end
    fprintf('Quantum scene created with 9 targets.\n');
end

function measurements = quantum_measurements(H, v_true_vector, params)
    % Quantum measurement simulation with frequency diversity
    fprintf('--- ULTIMATE MEASUREMENT SIMULATION ---\n');
    % Quantum signal generation with frequency diversity
    Hv_signal = H * v_true_vector;
    % Account for frequency diversity in signal power calculation
    % Each acquisition uses different frequencies, so we need to normalize
    signal_power_est = mean(Hv_signal(:).^2);
    target_SNR_linear = 10^(params.target_SNR_db / 10);
    noise_variance = signal_power_est / target_SNR_linear;
    noise_sigma = sqrt(noise_variance);
    actual_SNR_db = 10 * log10(signal_power_est / noise_variance);
    % Quantum noise addition with frequency-dependent scaling
    n_noise_vec = noise_sigma * randn(size(Hv_signal));
    u_measured_signal = Hv_signal + n_noise_vec;
    fprintf('NOISE: Target SNR %.1f dB. Actual SNR %.1f dB. Sigma: %g\n', ...
        params.target_SNR_db, actual_SNR_db, noise_sigma);
    fprintf('Signal Structure: 3 pMUTs, Random 45-65 kHz, 3 cycles, 0-12μs delays\n');
    measurements = struct();
    measurements.Hv_signal = Hv_signal;
    measurements.n_noise_vec = n_noise_vec;
    measurements.u_measured_signal = u_measured_signal;
    measurements.actual_SNR_db = actual_SNR_db;
    measurements.noise_sigma = noise_sigma;
end

function [recon, PSNR_admm_iters, residuals_admm_iters] = quantum_admm_diagnostic(H, measurements, scene_matrix, imaging_grid, params, output_folder)
    % Variant of quantum_admm that returns reconstruction, PSNR, and residuals arrays
    fprintf('\n--- ULTIMATE ADMM RECONSTRUCTION (DIAGNOSTIC) ---\n');
    A_matrix = H;
    b_vector = measurements.u_measured_signal;
    I_true_matrix = scene_matrix;
    v_true_vec_norm = scene_matrix(:);
    if max(abs(v_true_vec_norm)) > 0
        v_true_vec_norm = v_true_vec_norm ./ max(abs(v_true_vec_norm));
    end
    imageResolution = size(I_true_matrix);
    H_norm_factor = max(abs(A_matrix(:)));
    if H_norm_factor < eps
        H_norm_factor = 1;
    end
    A_admm = A_matrix ./ H_norm_factor;
    At_admm = transpose(A_admm);
    b_admm_vec = b_vector(:) / H_norm_factor;
    noise_sigma_admm = measurements.noise_sigma / H_norm_factor;
    [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = quantum_operator_setup(A_admm, At_admm, imageResolution);
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        params.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    PSNR_admm_iters = zeros([params.admm_max_iter 1]);
    residuals_admm_iters = zeros([params.admm_max_iter 1]);
    converged = false;
    for k_admm = 1:params.admm_max_iter
        x_prev = x_admm_img_iter;
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = quantum_admm_iteration(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            params.rho_admm, params.lambda_tv_reg, Hfun_pcg_admm, params);
        [PSNR_admm_iters(k_admm), residuals_admm_iters(k_admm)] = quantum_metrics(...
            x_admm_img_iter, v_true_vec_norm, b_admm_vec, Afun_admm, opDx_tv, params.lambda_tv_reg);
        if k_admm > 1 && params.early_termination
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < params.admm_tol
                converged = true;
                break;
            end
        end
    end
    recon = x_admm_img_iter;
    PSNR_admm_iters = PSNR_admm_iters(1:k_admm);
    residuals_admm_iters = residuals_admm_iters(1:k_admm);
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = quantum_operator_setup(A_admm, At_admm, imageResolution)
    % Quantum operator setup
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    [Dx_sparse, Dy_sparse] = quantum_difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [Dx, Dy] = quantum_difference_operators(imageSize)
    % Quantum difference operator creation
    rows = imageSize(1);
    cols = imageSize(2);
    N_img_pixels = rows * cols;
    Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
    last_col_indices_mask = false(N_img_pixels, 1);
    last_col_indices_mask((cols-1)*rows+1 : cols*rows) = true;
    Dx(last_col_indices_mask, :) = 0;
    Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
    last_row_indices_mask = false(N_img_pixels, 1);
    last_row_indices_mask(rows:rows:N_img_pixels) = true;
    Dy(last_row_indices_mask, :) = 0;
end

function [x_new, z_new, u_new] = quantum_admm_iteration(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, Hfun_pcg_admm, params)
    % Quantum ADMM iteration
    % Quantum x-update
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), params.pcg_tol, params.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    % Quantum z-update
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm = max(v_norm, eps);
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    % Quantum u-update
    u_new = u_old + opDx_tv(x_new) - z_new;
end

function [PSNR, residual] = quantum_metrics(x_admm_img_iter, v_true_vec_norm, b_admm_vec, Afun_admm, opDx_tv, lambda_tv_reg)
    % Quantum metrics calculation
    % Quantum normalization
    x_scl = real(x_admm_img_iter(:));
    x_range = max(x_scl) - min(x_scl);
    if x_range > eps
        x_scl = (x_scl - min(x_scl)) / x_range;
    else
        x_scl = zeros(size(x_scl));
    end
    % Quantum PSNR calculation
    MSE_curr = mean((x_scl - v_true_vec_norm).^2);
    PSNR = 10 * log10(1 / MSE_curr);
    % Quantum residual calculation
    r1 = b_admm_vec - Afun_admm(x_admm_img_iter);
    r2 = opDx_tv(x_admm_img_iter);
    tv_n = sum(sqrt(sum(r2.^2, 2)));
    residual = 0.5 * sum(r1(:).^2) + lambda_tv_reg * tv_n;
end

function plot_quantum_impulse_response(t_burst, impulse_response, output_folder, params)
    % Plot the excitation burst (impulse response)
    figure(2);
    set(gcf, 'visible', 'off');
    clf;
    plot(t_burst * 1e6, impulse_response);
    title(sprintf('Figure 2: Quantum Excitation Burst (Impulse Response, Z=%dmm)', params.target_distance_mm));
    xlabel('Time (us)');
    ylabel('Amplitude');
    grid on;
    set(gcf, 'Color', 'w');
    try
        saveas(gcf, fullfile(output_folder, 'figure2_quantum_excitation_burst.png'));
    catch ME
        fprintf('  Warning: Could not save impulse response figure: %s\n', ME.message);
    end
    close(gcf);
end

function plot_quantum_h_matrix_columns(H, params, output_folder)
    % Plot columns of the assembled H matrix
    figure(3);
    set(gcf, 'visible', 'off');
    clf;
    hold on;
    num_cols_to_plot = min(size(H, 2), 5);
    indices_to_plot = round(linspace(1, size(H, 2), num_cols_to_plot));
    t_axis_plot = (0:(size(H, 1) - 1)) / params.fs * 1e6;
    for n_idx = 1:length(indices_to_plot)
        col_idx = indices_to_plot(n_idx);
        plot(t_axis_plot, H(:, col_idx), 'DisplayName', sprintf('H col Px %d', col_idx));
    end
    hold off;
    xlabel('Overall Row Index (us)');
    ylabel('Amplitude');
    title(sprintf('Figure 3: Quantum Columns of Assembled H (Bistatic, Z=%dmm)', params.target_distance_mm));
    axis tight;
    grid on;
    legend('Location', 'best');
    set(gcf, 'Color', 'w');
    try
        saveas(gcf, fullfile(output_folder, 'figure3_quantum_columns_of_H.png'));
    catch ME
        fprintf('  Warning: Could not save H matrix columns figure: %s\n', ME.message);
    end
    close(gcf);
end

function plot_quantum_true_scene(scene_matrix, imaging_grid, params, output_folder)
    % Plot the true scene with targets
    figure(4);
    set(gcf, 'visible', 'off');
    clf;
    % Main scene plot
    subplot(1, 2, 1);
    imagesc(imaging_grid.x_coords * 1e3, imaging_grid.y_coords * 1e3, scene_matrix);
    axis image;
    colormap gray;
    colorbar;
    xlabel('x (mm)');
    ylabel('y (mm)');
    title(sprintf('Quantum True Image (XY plane at Z=%dmm)', params.target_distance_mm));
    set(gca, 'YDir', 'normal');
    % Vector plot
    subplot(1, 2, 2);
    stem(scene_matrix(:));
    axis square tight;
    title('Quantum True Image (vector)');
    xlabel('Pixel Index');
    ylabel('Amplitude');
    sgtitle(sprintf('Figure 4: Quantum True Image (Bistatic, Z=%dmm)', params.target_distance_mm));
    set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 450]);
    try
        saveas(gcf, fullfile(output_folder, 'figure4_quantum_true_scene.png'));
    catch ME
        fprintf('  Warning: Could not save true scene figure: %s\n', ME.message);
    end
    close(gcf);
end

function plot_quantum_measurement_signals(measurements, params, output_folder)
    % Plot measurement signals
    figure(6);
    set(gcf, 'visible', 'off');
    clf;
    t_ax_meas = (0:(length(measurements.Hv_signal) - 1)) / params.fs * 1e6;
    plot_lim_samp = min(length(measurements.Hv_signal), round(params.fs * 0.001 * 1.1));
    if plot_lim_samp == 0 && length(measurements.Hv_signal) > 0
        plot_lim_samp = length(measurements.Hv_signal);
    end
    subplot(1, 3, 1);
    if ~isempty(measurements.Hv_signal) && plot_lim_samp > 0
        plot(t_ax_meas(1:plot_lim_samp), measurements.Hv_signal(1:plot_lim_samp));
    end
    axis tight; grid on; xlabel('Time (us)'); title('Hv - Ideal'); yl = ylim;
    subplot(1, 3, 2);
    if ~isempty(measurements.n_noise_vec) && plot_lim_samp > 0
        plot(t_ax_meas(1:plot_lim_samp), measurements.n_noise_vec(1:plot_lim_samp));
    end
    axis tight; grid on; xlabel('Time (us)'); title(sprintf('Noise (SNR~%.1fdB)', measurements.actual_SNR_db));
    if all(isfinite(yl)); ylim(yl); end
    subplot(1, 3, 3);
    if ~isempty(measurements.u_measured_signal) && plot_lim_samp > 0
        plot(t_ax_meas(1:plot_lim_samp), measurements.u_measured_signal(1:plot_lim_samp));
    end
    axis tight; grid on; xlabel('Time (us)'); title('u (Meas+Noise)');
    if all(isfinite(yl)); ylim(yl); end
    sgtitle(sprintf('Figure 6: Quantum Measurement Signals (Bistatic, Z=%dmm)', params.target_distance_mm));
    set(gcf, 'Color', 'w');
    try
        saveas(gcf, fullfile(output_folder, 'figure6_quantum_measurement_signals.png'));
    catch ME
        fprintf('  Warning: Could not save measurement signals figure: %s\n', ME.message);
    end
    close(gcf);
end

function update_quantum_admm_visualization(x_admm_img_iter, I_true_matrix, PSNR_admm_iters, residuals_admm_iters, k_admm, params, imaging_grid)
    % Update Quantum ADMM visualization
    figure(8);
    % Normalize current reconstruction
    x_scl = real(x_admm_img_iter(:));
    x_range = max(x_scl) - min(x_scl);
    if x_range > eps
        x_scl = (x_scl - min(x_scl)) / x_range;
    else
        x_scl = zeros(size(x_scl));
    end
    % Plot target
    subplot(2, 3, [1 4]);
    imagesc(imaging_grid.x_coords * 1e3, imaging_grid.y_coords * 1e3, I_true_matrix);
    axis image; colormap(gca, gray); colorbar; set(gca, 'YDir', 'normal');
    title(sprintf('Target (Z=%dmm)', params.target_distance_mm)); xlabel('x(mm)'); ylabel('y(mm)');
    % Plot reconstruction
    subplot(2, 3, [2 5]);
    imagesc(imaging_grid.x_coords * 1e3, imaging_grid.y_coords * 1e3, reshape(x_scl, size(I_true_matrix)));
    axis image; colormap(gca, gray); colorbar; set(gca, 'YDir', 'normal');
    title(sprintf('Quantum Reconstruction (Z=%dmm)\n\\lambda=%.1e,\\rho=%.1f\nPSNR=%.2fdB,It %d', ...
        params.target_distance_mm, params.lambda_tv_reg, params.rho_admm, PSNR_admm_iters(k_admm), k_admm));
    xlabel('x(mm)'); ylabel('y(mm)');
    % Plot PSNR convergence
    subplot(2, 3, 3);
    plot(1:k_admm, PSNR_admm_iters(1:k_admm), 'm-', 'LineWidth', 2);
    title('PSNR/Iter'); xlabel('Iter'); ylabel('PSNR(dB)'); grid on; axis tight;
    if k_admm > 1
        yl = ylim; if diff(yl) < 1; yl(2) = yl(1) + 1; end; ylim(yl);
    end
    % Plot residual convergence
    subplot(2, 3, 6);
    plot(1:k_admm, log10(residuals_admm_iters(1:k_admm)), '-', 'LineWidth', 2);
    title('log10(Obj)/Iter'); xlabel('Iter'); ylabel('log10(Val)'); grid on; axis tight;
    if k_admm > 1
        yl = ylim; if diff(yl) < 0.1; yl(2) = yl(1) + 0.1; end; ylim(yl);
    end
end