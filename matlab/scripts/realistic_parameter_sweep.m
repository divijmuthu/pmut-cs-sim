%% realistic_parameter_sweep.m - Comprehensive Parameter Sweep for Realistic pMUT Reconstruction
% Sweeps through various parameters and ranks results by PSNR
% Tests grid resolution, target sizes, spacing, acquisitions, and ADMM parameters

clear; clc; close all;

%% ===== SWEEP CONFIGURATION =====
fprintf('=== REALISTIC pMUT PARAMETER SWEEP ===\n');
fprintf('Testing various configurations and ranking by PSNR\n\n');

% Base configuration (best from V29 sweep)
base_config = struct();
base_config.c = 343;                    % Speed of sound (m/s)
base_config.fs = 1e6;                   % Sampling frequency (Hz)
base_config.pmut_width_m = 0.020;       % pMUT width (m)
base_config.tx_pool_width_m = 0.200;    % Transmitter pool width (m)
base_config.grid_width_m = 0.150;       % Imaging grid width (m)
base_config.target_distance_m = 0.150;  % Target distance (m)
base_config.grid_depth_range_m = 0.100; % Grid depth range (m)
base_config.excitation_amplitude = 1e15; % Excitation amplitude

% REALISTIC pMUT PARAMETERS (from experimental data)
base_config.pmut_resonance_freq = 57700; % 57.7 kHz average resonance
base_config.pmut_bandwidth = 2520;      % 2.52 kHz standard deviation
base_config.impulse_duration_us = 10;   % Short impulse excitation (10 μs)

% BEST PERFORMING PARAMETERS from V29 sweep
base_config.num_active_tx = 5;          % Best from V28_final
base_config.max_delay_rand_us = 500;    % Best from V28_final
base_config.apodization_mode = 'uniform'; % Best from V28_final
base_config.frequency_offset_hz = 0;    % No frequency offset

% ADMM Reconstruction Parameters (from TheoryTest2)
base_config.numItersADMM = 50;          % ADMM iterations
base_config.rho_admm = 6.73;            % Optimized ADMM penalty
base_config.lambda_tv_reg = 0.7438;     % Optimized TV regularization
base_config.admm_tol = 1.2e-5;          % Optimized tolerance
base_config.admm_max_iter = 50;         % Fixed max iterations
base_config.pcg_max_iter = 30;          % Reduced PCG iterations
base_config.pcg_tol = 1e-8;             % Slightly relaxed PCG tolerance
base_config.target_SNR_db = 35;         % Target SNR for noise

%% ===== PARAMETER SWEEP DEFINITIONS =====
sweep_params = struct();

% Grid resolution sweep (balance resolution vs computation)
sweep_params.grid_step_mm = [5, 7.5, 10];  % Grid step in mm

% Target size sweep (test spatial resolution)
sweep_params.target_size_mm = [3, 4, 6, 8];  % Target diameter in mm

% Grid spacing sweep (test target separation)
sweep_params.grid_spacing_mm = [15, 20, 25];  % Spacing between targets

% Number of acquisitions sweep (test measurement diversity)
sweep_params.num_acquisitions = [15, 20, 25, 30];  % Number of acquisitions

% ADMM parameter sweep (test reconstruction quality)
sweep_params.lambda_tv_reg = [0.5, 0.7438, 1.0, 1.5];  % TV regularization strength

% Target pattern sweep (test different arrangements)
sweep_params.target_patterns = {'3x3_grid', '2x2_grid', 'line_5', 'cross_5'};

%% ===== OUTPUT SETUP =====
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('realistic_parameter_sweep_output', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('Output folder: %s\n\n', output_folder);

%% ===== GENERATE PARAMETER COMBINATIONS =====
% Create parameter combinations (focus on most important ones first)
param_combinations = [];

% Primary sweep: Grid resolution + Target size + Grid spacing
for grid_step = sweep_params.grid_step_mm
    for target_size = sweep_params.target_size_mm
        for grid_spacing = sweep_params.grid_spacing_mm
            % Skip combinations that don't make sense
            if target_size >= grid_spacing
                continue; % Target too big for spacing
            end
            
            param_combinations = [param_combinations; struct(...
                'grid_step_mm', grid_step, ...
                'target_size_mm', target_size, ...
                'grid_spacing_mm', grid_spacing, ...
                'num_acquisitions', 20, ... % Fixed for primary sweep
                'lambda_tv_reg', 0.7438, ... % Fixed for primary sweep
                'target_pattern', '3x3_grid' ... % Fixed for primary sweep
            )];
        end
    end
end

% Secondary sweep: Acquisitions + ADMM parameters
for num_acq = sweep_params.num_acquisitions
    for lambda_tv = sweep_params.lambda_tv_reg
        param_combinations = [param_combinations; struct(...
            'grid_step_mm', 7.5, ... % Fixed for secondary sweep
            'target_size_mm', 4, ... % Fixed for secondary sweep
            'grid_spacing_mm', 20, ... % Fixed for secondary sweep
            'num_acquisitions', num_acq, ...
            'lambda_tv_reg', lambda_tv, ...
            'target_pattern', '3x3_grid' ...
        )];
    end
end

% Pattern sweep: Different target arrangements
for pattern = sweep_params.target_patterns
    param_combinations = [param_combinations; struct(...
        'grid_step_mm', 7.5, ... % Fixed for pattern sweep
        'target_size_mm', 4, ... % Fixed for pattern sweep
        'grid_spacing_mm', 20, ... % Fixed for pattern sweep
        'num_acquisitions', 20, ... % Fixed for pattern sweep
        'lambda_tv_reg', 0.7438, ... % Fixed for pattern sweep
        'target_pattern', pattern{1} ...
    )];
end

num_combinations = length(param_combinations);
fprintf('Total parameter combinations: %d\n', num_combinations);
fprintf('Estimated runtime: ~%d minutes\n\n', round(num_combinations * 2));

%% ===== RESULTS STORAGE =====
results = struct();
results.configs = cell(num_combinations, 1);
results.metrics = cell(num_combinations, 1);
results.times = cell(num_combinations, 1);
results.psnr_rankings = zeros(num_combinations, 1);

%% ===== MAIN SWEEP LOOP =====
for i = 1:num_combinations
    fprintf('\n=== [%d/%d] PARAMETER COMBINATION ===\n', i, num_combinations);
    
    % Extract parameters
    params = param_combinations(i);
    fprintf('Grid step: %.1f mm, Target size: %.1f mm, Grid spacing: %.1f mm\n', ...
        params.grid_step_mm, params.target_size_mm, params.grid_spacing_mm);
    fprintf('Acquisitions: %d, Lambda TV: %.4f, Pattern: %s\n', ...
        params.num_acquisitions, params.lambda_tv_reg, params.target_pattern);
    
    % Create configuration
    config = base_config;
    config.grid_step_m = params.grid_step_mm / 1000;
    config.num_acquisitions = params.num_acquisitions;
    config.lambda_tv_reg = params.lambda_tv_reg;
    
    % Create run name
    run_name = sprintf('run%03d_gs%.1f_ts%.1f_sp%.1f_acq%d_ltv%.4f_%s', ...
        i, params.grid_step_mm, params.target_size_mm, params.grid_spacing_mm, ...
        params.num_acquisitions, params.lambda_tv_reg, params.target_pattern);
    
    fprintf('Run name: %s\n', run_name);
    
    % Run reconstruction
    try
        tic;
        [psnr, correlation, reconstruction_time, scene_info] = run_single_reconstruction(config, params, run_name, output_folder);
        total_time = toc;
        
        % Store results
        results.configs{i} = params;
        results.metrics{i} = struct('psnr', psnr, 'correlation', correlation, 'scene_info', scene_info);
        results.times{i} = struct('reconstruction_time', reconstruction_time, 'total_time', total_time);
        results.psnr_rankings(i) = psnr;
        
        fprintf('Results: PSNR=%.2f dB, Correlation=%.4f, Time=%.1f s\n', ...
            psnr, correlation, total_time);
        
    catch ME
        fprintf('ERROR in run %d: %s\n', i, ME.message);
        results.configs{i} = params;
        results.metrics{i} = struct('psnr', -inf, 'correlation', 0, 'scene_info', struct());
        results.times{i} = struct('reconstruction_time', 0, 'total_time', 0);
        results.psnr_rankings(i) = -inf;
    end
end

%% ===== RANK AND ANALYZE RESULTS =====
fprintf('\n=== RANKING RESULTS BY PSNR ===\n');

% Sort by PSNR (descending)
[sorted_psnr, sort_indices] = sort(results.psnr_rankings, 'descend');

% Display top 10 results
fprintf('\nTop 10 Results by PSNR:\n');
fprintf('Rank | PSNR (dB) | Correlation | Grid Step | Target Size | Grid Spacing | Acquisitions | Lambda TV | Pattern | Time (s)\n');
fprintf('-----|-----------|------------|----------|-------------|--------------|-------------|-----------|---------|---------\n');

for rank = 1:min(10, num_combinations)
    idx = sort_indices(rank);
    params = results.configs{idx};
    metrics = results.metrics{idx};
    times = results.times{idx};
    
    fprintf('%4d | %9.2f | %10.4f | %8.1f | %11.1f | %13.1f | %12d | %9.4f | %7s | %8.1f\n', ...
        rank, metrics.psnr, metrics.correlation, params.grid_step_mm, params.target_size_mm, ...
        params.grid_spacing_mm, params.num_acquisitions, params.lambda_tv_reg, ...
        params.target_pattern, times.total_time);
end

% Display worst 5 results
fprintf('\nWorst 5 Results by PSNR:\n');
fprintf('Rank | PSNR (dB) | Correlation | Grid Step | Target Size | Grid Spacing | Acquisitions | Lambda TV | Pattern | Time (s)\n');
fprintf('-----|-----------|------------|----------|-------------|--------------|-------------|-----------|---------|---------\n');

for rank = max(1, num_combinations-4):num_combinations
    idx = sort_indices(rank);
    params = results.configs{idx};
    metrics = results.metrics{idx};
    times = results.times{idx};
    
    fprintf('%4d | %9.2f | %10.4f | %8.1f | %11.1f | %13.1f | %12d | %9.4f | %7s | %8.1f\n', ...
        rank, metrics.psnr, metrics.correlation, params.grid_step_mm, params.target_size_mm, ...
        params.grid_spacing_mm, params.num_acquisitions, params.lambda_tv_reg, ...
        params.target_pattern, times.total_time);
end

%% ===== SAVE COMPLETE RESULTS =====
% Save all results
save(fullfile(output_folder, 'complete_sweep_results.mat'), 'results', 'base_config', 'sweep_params', 'param_combinations');

% Create summary CSV
create_summary_csv(results, sort_indices, output_folder);

% Create analysis plots
create_analysis_plots(results, sort_indices, output_folder);

fprintf('\n=== PARAMETER SWEEP COMPLETE ===\n');
fprintf('Results saved to: %s\n', output_folder);
fprintf('Best PSNR: %.2f dB\n', max(results.psnr_rankings));
fprintf('Worst PSNR: %.2f dB\n', min(results.psnr_rankings));
fprintf('Average PSNR: %.2f dB\n', mean(results.psnr_rankings));

%% ===== HELPER FUNCTIONS =====

function [psnr, correlation, reconstruction_time, scene_info] = run_single_reconstruction(config, params, run_name, output_folder)
    % Run a single reconstruction with given parameters
    
    % Create subfolder for this run
    run_folder = fullfile(output_folder, run_name);
    if ~exist(run_folder, 'dir')
        mkdir(run_folder);
    end
    
    % Generate H matrix
    fprintf('  Generating H matrix...\n');
    [H, coherence_matrix] = generate_h_matrix_realistic_pmut(config);
    
    % Create imaging grid
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    
    % Create scene based on pattern
    scene_matrix = create_scene_with_pattern(X_mesh, Z_mesh, config, params);
    v_true_vector = scene_matrix(:);
    
    % Forward simulation
    fprintf('  Running forward simulation...\n');
    Hv_signal = H * v_true_vector;
    signal_power = mean(Hv_signal(:).^2);
    target_SNR_linear = 10^(config.target_SNR_db / 10);
    noise_variance = signal_power / target_SNR_linear;
    noise_sigma = sqrt(noise_variance);
    noise = noise_sigma * randn(size(Hv_signal));
    b_vector = Hv_signal + noise;
    
    % ADMM reconstruction
    fprintf('  Running ADMM reconstruction...\n');
    tic;
    reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, run_folder);
    reconstruction_time = toc;
    
    % Calculate metrics
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    scene_vec = scene_norm(:);
    recon_vec = recon_norm(:);
    correlation = (scene_vec - mean(scene_vec))' * (recon_vec - mean(recon_vec)) / ...
        (sqrt(sum((scene_vec - mean(scene_vec)).^2)) * sqrt(sum((recon_vec - mean(recon_vec)).^2)));
    
    % Scene info
    scene_info = struct();
    scene_info.grid_size = size(scene_matrix);
    scene_info.num_pixels = numel(scene_matrix);
    scene_info.num_targets = nnz(scene_matrix);
    scene_info.h_matrix_size = size(H);
    
    % Save run results
    save(fullfile(run_folder, 'run_results.mat'), 'psnr', 'correlation', 'reconstruction_time', 'scene_info');
end

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

function create_summary_csv(results, sort_indices, output_folder)
    % Create summary CSV file
    csv_data = [];
    
    for i = 1:length(results.configs)
        idx = sort_indices(i);
        params = results.configs{idx};
        metrics = results.metrics{idx};
        times = results.times{idx};
        
        row = struct();
        row.rank = i;
        row.psnr = metrics.psnr;
        row.correlation = metrics.correlation;
        row.grid_step_mm = params.grid_step_mm;
        row.target_size_mm = params.target_size_mm;
        row.grid_spacing_mm = params.grid_spacing_mm;
        row.num_acquisitions = params.num_acquisitions;
        row.lambda_tv_reg = params.lambda_tv_reg;
        row.target_pattern = params.target_pattern;
        row.reconstruction_time = times.reconstruction_time;
        row.total_time = times.total_time;
        
        csv_data = [csv_data; row];
    end
    
    csv_table = struct2table(csv_data);
    writetable(csv_table, fullfile(output_folder, 'parameter_sweep_summary.csv'));
end

function create_analysis_plots(results, sort_indices, output_folder)
    % Create analysis plots
    
    % Extract data
    psnr_values = zeros(length(results.configs), 1);
    grid_steps = zeros(length(results.configs), 1);
    target_sizes = zeros(length(results.configs), 1);
    grid_spacings = zeros(length(results.configs), 1);
    num_acquisitions = zeros(length(results.configs), 1);
    lambda_tv_values = zeros(length(results.configs), 1);
    
    for i = 1:length(results.configs)
        psnr_values(i) = results.metrics{i}.psnr;
        grid_steps(i) = results.configs{i}.grid_step_mm;
        target_sizes(i) = results.configs{i}.target_size_mm;
        grid_spacings(i) = results.configs{i}.grid_spacing_mm;
        num_acquisitions(i) = results.configs{i}.num_acquisitions;
        lambda_tv_values(i) = results.configs{i}.lambda_tv_reg;
    end
    
    % PSNR vs parameters plots
    figure('Visible', 'off');
    
    subplot(2, 3, 1);
    scatter(grid_steps, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Step (mm)');
    ylabel('PSNR (dB)');
    title('PSNR vs Grid Step');
    grid on;
    
    subplot(2, 3, 2);
    scatter(target_sizes, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Target Size (mm)');
    ylabel('PSNR (dB)');
    title('PSNR vs Target Size');
    grid on;
    
    subplot(2, 3, 3);
    scatter(grid_spacings, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Spacing (mm)');
    ylabel('PSNR (dB)');
    title('PSNR vs Grid Spacing');
    grid on;
    
    subplot(2, 3, 4);
    scatter(num_acquisitions, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Number of Acquisitions');
    ylabel('PSNR (dB)');
    title('PSNR vs Acquisitions');
    grid on;
    
    subplot(2, 3, 5);
    scatter(lambda_tv_values, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Lambda TV');
    ylabel('PSNR (dB)');
    title('PSNR vs Lambda TV');
    grid on;
    
    subplot(2, 3, 6);
    histogram(psnr_values, 20);
    xlabel('PSNR (dB)');
    ylabel('Count');
    title('PSNR Distribution');
    grid on;
    
    saveas(gcf, fullfile(output_folder, 'parameter_analysis.png'));
    close(gcf);
end

% Include all the helper functions from realistic_reconstruction_demo.m
% (generate_h_matrix_realistic_pmut, run_admm_reconstruction, etc.)

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