%% sweepScriptV25.m - Fixed H-Matrix Optimization with Time-Domain Diversity
% Combines V24's advanced analytics with V23's working time-domain approach
% Based on successful results from V23 showing 0.791 max coherence

clear; clc; close all;

% Add paths
addpath(genpath(pwd));
addpath(genpath(fullfile(pwd, 'm_files')));

%% Base Configuration (Optimized from V23 results)
base_config = struct();
base_config.c = 343;                    % Speed of sound (m/s)
base_config.fs = 1e6;                   % Sampling frequency (Hz)
base_config.pmut_width_m = 0.020;       % pMUT width (m)
base_config.tx_pool_width_m = 0.200;    % Transmitter pool width (m)
base_config.grid_width_m = 0.150;       % Imaging grid width (m)
base_config.target_distance_m = 0.150;  % Target distance (m)
base_config.grid_depth_range_m = 0.100; % Grid depth range (m)
base_config.grid_step_m = 0.010;        % Grid step (m)
base_config.num_acquisitions = 30;      % Increased for better statistics
base_config.excitation_amplitude = 1e15; % Excitation amplitude
base_config.f_min_hz = 45000;           % Min frequency (Hz)
base_config.f_max_hz = 65000;           % Max frequency (Hz)

%% Targeted Parameter Sweep (Based on V23 Best Results)
p = struct();

% Focus on best performers: 5 transmitters, 500μs delays, uniform apodization
p.num_active_tx_sweep = {5};            % Best from V23
p.delay_rand_sweep_us = {500};          % Best from V23
p.apodization_sweep = {'uniform'};      % Best from V23

% New targeted parameters to optimize further
p.grid_step_sweep = {0.008, 0.010, 0.012};  % Grid density optimization
p.tx_pool_width_sweep = {0.200, 0.250, 0.300}; % Transmitter spacing (avoiding kerf issues)
p.frequency_bandwidth_sweep = {[45000, 65000], [40000, 70000], [50000, 60000]}; % Frequency diversity
p.excitation_pattern_sweep = {'chirp', 'gaussian', 'hanning'}; % Excitation patterns

%% Advanced Analytics Configuration
analytics_config = struct();
analytics_config.compute_rip_proxy = true;      % RIP proxy calculation
analytics_config.compute_condition_number = true; % Condition number analysis
analytics_config.compute_svd_decay = true;      % SVD decay analysis
analytics_config.compute_correlation_matrix = true; % Correlation analysis
analytics_config.compute_energy_distribution = true; % Energy distribution
analytics_config.compute_rank_analysis = true;  % Rank analysis
analytics_config.save_detailed_plots = true;    % Save detailed plots
analytics_config.compute_compressed_sensing_metrics = true; % CS-specific metrics

%% Output Configuration
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('sweep_output_fixed_v25', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('--- Fixed H-Matrix Optimization V25 (with Time-Domain Diversity) ---\n');
fprintf('Saving results to: %s\n\n', output_folder);

%% Generate Parameter Combinations
param_grid = allcomb(p.num_active_tx_sweep, p.delay_rand_sweep_us, p.apodization_sweep, ...
                     p.grid_step_sweep, p.tx_pool_width_sweep, p.frequency_bandwidth_sweep, p.excitation_pattern_sweep);

num_combinations = size(param_grid, 1);
fprintf('Total parameter combinations: %d\n\n', num_combinations);

%% Results Storage
results = struct();
results.configs = cell(num_combinations, 1);
results.metrics = cell(num_combinations, 1);
results.H_matrices = cell(num_combinations, 1);
results.scenes = cell(num_combinations, 1);

%% Main Sweep Loop
for i = 1:num_combinations
    run_name = sprintf('run%03d_tx%d_del%d_%s_gs%.3f_tpw%.3f_fb%d_%s', i, ...
        param_grid{i, 1}, param_grid{i, 2}, param_grid{i, 3}, ...
        param_grid{i, 4}, param_grid{i, 5}, param_grid{i, 6}(2) - param_grid{i, 6}(1), param_grid{i, 7});
    fprintf('--- [%d/%d] Running Test: %s ---\n', i, num_combinations, run_name);
    
    % Extract parameters
    config = base_config;
    config.num_active_tx = param_grid{i, 1};
    config.max_delay_rand_us = param_grid{i, 2};
    config.apodization_mode = param_grid{i, 3};
    config.grid_step_m = param_grid{i, 4};
    config.tx_pool_width_m = param_grid{i, 5};
    config.f_min_hz = param_grid{i, 6}(1);
    config.f_max_hz = param_grid{i, 6}(2);
    config.excitation_pattern = param_grid{i, 7};
    
    % Generate H matrix using V23's time-domain approach
    [H, scene, acquisition_info] = generate_h_matrix_v23_style(config);
    
    % Run comprehensive diagnostics
    metrics = run_advanced_diagnostics(H, scene, acquisition_info, analytics_config);
    
    % Store results
    results.configs{i} = config;
    results.metrics{i} = metrics;
    results.H_matrices{i} = H;
    results.scenes{i} = scene;
    
    % Save individual run results
    run_folder = fullfile(output_folder, run_name);
    if ~exist(run_folder, 'dir')
        mkdir(run_folder);
    end
    
    % Save detailed results
    save(fullfile(run_folder, 'results.mat'), 'H', 'scene', 'config', 'metrics', 'acquisition_info');
    
    % Generate and save plots
    if analytics_config.save_detailed_plots
        generate_detailed_plots(H, scene, metrics, run_folder, run_name);
    end
    
    % Print summary
    fprintf('Max Coherence: %.4f, Sparsity: %.2f%%, Non-zero Columns: %d/%d\n', ...
        full(metrics.max_coherence), full(metrics.sparsity), full(metrics.num_nonzero_cols), full(metrics.total_cols));
    fprintf('Condition Number: %.2f, RIP Proxy: %.4f\n', ...
        full(metrics.condition_number), full(metrics.rip_proxy));
    fprintf('Energy Concentration: %.2f%%, Rank: %d/%d\n\n', ...
        full(metrics.energy_concentration), full(metrics.rank), min(size(H)));
end

%% Generate Summary Report
generate_summary_report(results, output_folder);

%% Save Complete Results
save(fullfile(output_folder, 'complete_results.mat'), 'results', 'base_config', 'p', 'analytics_config');

fprintf('--- Fixed Optimization Complete! ---\n');
fprintf('Results saved to: %s\n', output_folder);

%% Helper Functions

function [H, scene, acquisition_info] = generate_h_matrix_v23_style(config)
    % Initialize Field II
    field_init(-1);
    
    % Create scene
    scene = create_test_scene(config);
    
    % Setup Field II transducers (using 2D array approach like V23)
    fs = config.fs;
    c = config.c;
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
    
    % Create imaging grid
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    
    % Initialize acquisition storage
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    fprintf('Starting %d acquisitions...\n', config.num_acquisitions);
    tic;
    
    for r_acq = 1:config.num_acquisitions
        acq_tic = tic;
        
        % Randomly choose 3 or 5 transmitters for this acquisition
        if rand < 0.5
            num_active_tx = 3;
        else
            num_active_tx = 5;
        end
        
        % Generate diverse excitation signals
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        
        % Generate excitation pattern based on config
        excitation_pattern = generate_excitation_pattern_v25(config, num_active_tx);
        
        % Setup apodization
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random')
            apod_weights = rand(1, num_active_tx);
        end
        
        % Apply excitation pattern to Field II
        xdc_impulse(TxAperture, excitation_pattern);
        
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
        
        if mod(r_acq, 5) == 0 || r_acq == config.num_acquisitions
            fprintf('Acquisition %d/%d completed in %.2f seconds.\n', r_acq, config.num_acquisitions, toc(acq_tic));
        end
    end
    
    total_time = toc;
    fprintf('All acquisitions completed in %.2f seconds.\n', total_time);
    fprintf('Assembling H-matrix using interpolation...\n');
    
    % Assemble H matrix using V23's time-domain approach
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(config.num_acquisitions, N_pixels);
        acquisition_info = struct();
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
        fprintf('H matrix stats: min=%.3e, max=%.3e, nnz=%d, sum=%.3e\n', ...
            min(nonzeros(H)), max(nonzeros(H)), nnz(H), sum(nonzeros(H)));
    else
        fprintf('H matrix stats: all zeros, nnz=0\n');
    end
    
    % Column analysis
    col_norms = vecnorm(H, 2, 1);
    num_nz_cols = sum(col_norms > 1e-20);
    fprintf('Number of nonzero columns: %d (out of %d)\n', full(num_nz_cols), full(size(H,2)));
    
    % Store acquisition info
    acquisition_info = struct();
    acquisition_info.total_rows = total_rows;
    acquisition_info.K_global_per_acq = K_global_per_acq;
    acquisition_info.t_common_axis = t_common_axis;
end

function scene = create_test_scene(config)
    % Create imaging grid
    x_coords = -config.grid_width_m/2:config.grid_step_m:config.grid_width_m/2;
    z_coords = config.target_distance_m - config.grid_depth_range_m/2:config.grid_step_m:config.target_distance_m + config.grid_depth_range_m/2;
    [X, Z] = meshgrid(x_coords, z_coords);
    
    scene.x_coords = X(:);
    scene.z_coords = Z(:);
    scene.N_pixels = length(scene.x_coords);
    
    % Create test scene (point targets)
    scene.targets = struct();
    scene.targets.x = [0, 0.02, -0.02, 0.01, -0.01];
    scene.targets.z = [config.target_distance_m, config.target_distance_m + 0.02, config.target_distance_m - 0.02, config.target_distance_m + 0.01, config.target_distance_m - 0.01];
    scene.targets.amplitude = [1.0, 0.8, 0.8, 0.6, 0.6];
end

function excitation_pattern = generate_excitation_pattern_v25(config, num_active_tx)
    % Generate excitation pattern based on configuration
    switch config.excitation_pattern
        case 'chirp'
            % Chirp excitation
            t = linspace(0, 1e-6, 100);
            f0 = config.f_min_hz;
            f1 = config.f_max_hz;
            excitation_pattern = chirp(t, f0, 1e-6, f1);
        case 'gaussian'
            % Gaussian pulse
            t = linspace(-2e-6, 2e-6, 100);
            sigma = 0.5e-6;
            excitation_pattern = exp(-(t.^2)/(2*sigma^2));
        case 'hanning'
            % Hanning window
            t = linspace(0, 1e-6, 100);
            excitation_pattern = 0.5 * (1 - cos(2*pi*t/1e-6));
        otherwise
            % Default: simple pulse
            excitation_pattern = ones(1, 50);
    end
    
    % Normalize and apply excitation amplitude
    excitation_pattern = excitation_pattern / max(abs(excitation_pattern));
    excitation_pattern = excitation_pattern * config.excitation_amplitude * (0.8 + 0.4*rand());
end

function h_aligned = interpolate_h_response(t_current, h_current, t_common)
    % Interpolate H response to common time axis (from V23)
    N_pixels = size(h_current, 2);
    h_aligned = zeros(length(t_common), N_pixels);
    if length(t_current) > 1
        for px_col = 1:N_pixels
            h_aligned(:, px_col) = interp1(t_current, h_current(:, px_col), t_common, 'linear', 0);
        end
    end
end

function metrics = run_advanced_diagnostics(H, scene, acquisition_info, analytics_config)
    % Basic metrics
    metrics = struct();
    
    % Sparsity
    metrics.sparsity = (1 - nnz(H) / numel(H)) * 100;
    
    % Column analysis
    col_norms = vecnorm(H, 2, 1);
    metrics.num_nonzero_cols = sum(col_norms > 1e-20);
    metrics.total_cols = size(H, 2);
    metrics.column_coverage = metrics.num_nonzero_cols / metrics.total_cols;
    
    % Coherence analysis
    H_normalized = H ./ vecnorm(H, 2, 1);
    coherence_matrix = abs(H_normalized' * H_normalized);
    coherence_matrix(eye(size(coherence_matrix)) == 1) = 0; % Remove diagonal
    metrics.max_coherence = max(coherence_matrix(:));
    metrics.mean_coherence = mean(coherence_matrix(:));
    
    % Advanced metrics
    if analytics_config.compute_condition_number
        metrics.condition_number = cond(full(H));
    end
    
    if analytics_config.compute_rip_proxy
        % RIP proxy using random submatrices
        K = min(10, size(H, 2));
        rip_values = zeros(5, 1);
        for i = 1:5
            idx = randperm(size(H, 2), K);
            s = svd(full(H(:, idx)));
            rip_values(i) = max(s) / min(s);
        end
        metrics.rip_proxy = mean(rip_values);
    end
    
    if analytics_config.compute_rank_analysis
        metrics.rank = rank(full(H));
        metrics.rank_ratio = metrics.rank / min(size(H));
    end
    
    if analytics_config.compute_energy_distribution
        % Energy concentration in top singular values
        [~, S, ~] = svd(full(H), 'econ');
        S = diag(S);
        total_energy = sum(S.^2);
        top_energy = sum(S(1:min(10, length(S))).^2);
        metrics.energy_concentration = (top_energy / total_energy) * 100;
    end
    
    % Compressed sensing metrics
    if analytics_config.compute_compressed_sensing_metrics
        metrics.spark_lower_bound = estimate_spark_lower_bound(H);
        metrics.nullspace_property = check_nullspace_property(H);
    end
end

function spark_lb = estimate_spark_lower_bound(H)
    % Estimate lower bound on spark (minimum number of linearly dependent columns)
    [~, S, ~] = svd(full(H), 'econ');
    S = diag(S);
    spark_lb = sum(S > 1e-12) + 1;
end

function nullspace_prop = check_nullspace_property(H)
    % Check if H satisfies nullspace property for sparse recovery
    % This is a simplified check
    nullspace_prop = rank(full(H)) == size(H, 1);
end

function generate_detailed_plots(H, scene, metrics, output_folder, run_name)
    % SVD decay plot
    [~, S, ~] = svd(full(H), 'econ');
    S = diag(S);
    
    fig_svd = figure('Name', ['SVD Decay: ' run_name], 'visible', 'off');
    semilogy(S, 'LineWidth', 2);
    title('Singular Value Decay');
    xlabel('Singular Value Index'); ylabel('Magnitude');
    grid on;
    saveas(fig_svd, fullfile(output_folder, ['svd_decay_' run_name '.png']));
    close(fig_svd);
    
    % Coherence histogram
    H_normalized = H ./ vecnorm(H, 2, 1);
    coherence_matrix = abs(H_normalized' * H_normalized);
    coherence_matrix(eye(size(coherence_matrix)) == 1) = 0;
    coherence_values = coherence_matrix(coherence_matrix > 0);
    
    fig_coherence = figure('Name', ['Coherence Histogram: ' run_name], 'visible', 'off');
    histogram(coherence_values, 50);
    title('Mutual Coherence Distribution');
    xlabel('Coherence Value'); ylabel('Frequency');
    grid on;
    saveas(fig_coherence, fullfile(output_folder, ['coherence_hist_' run_name '.png']));
    close(fig_coherence);
    
    % H matrix visualization
    fig_H = figure('Name', ['H Matrix: ' run_name], 'visible', 'off');
    imagesc(full(H));
    colorbar;
    title('H Matrix Visualization');
    xlabel('Pixel Index'); ylabel('Acquisition Index');
    saveas(fig_H, fullfile(output_folder, ['H_matrix_' run_name '.png']));
    close(fig_H);
end

function generate_summary_report(results, output_folder)
    % Create summary CSV
    summary_data = [];
    
    for i = 1:length(results.metrics)
        metrics = results.metrics{i};
        config = results.configs{i};
        
        row = struct();
        row.run_index = i;
        row.max_coherence = metrics.max_coherence;
        row.mean_coherence = metrics.mean_coherence;
        row.sparsity = metrics.sparsity;
        row.column_coverage = metrics.column_coverage;
        row.condition_number = metrics.condition_number;
        row.rip_proxy = metrics.rip_proxy;
        row.rank_ratio = metrics.rank_ratio;
        row.energy_concentration = metrics.energy_concentration;
        
        % Configuration parameters
        row.num_active_tx = config.num_active_tx;
        row.max_delay_rand_us = config.max_delay_rand_us;
        row.grid_step_m = config.grid_step_m;
        row.tx_pool_width_m = config.tx_pool_width_m;
        row.f_bandwidth = config.f_max_hz - config.f_min_hz;
        row.excitation_pattern = config.excitation_pattern;
        
        summary_data = [summary_data; row];
    end
    
    % Convert to table and save
    summary_table = struct2table(summary_data);
    writetable(summary_table, fullfile(output_folder, 'detailed_summary.csv'));
    
    % Print best performers
    fprintf('\n--- Best Performers ---\n');
    coherence_values = zeros(length(results.metrics), 1);
    coverage_values = zeros(length(results.metrics), 1);
    for i = 1:length(results.metrics)
        coherence_values(i) = results.metrics{i}.max_coherence;
        coverage_values(i) = results.metrics{i}.column_coverage;
    end
    [~, best_coherence_idx] = min(coherence_values);
    [~, best_coverage_idx] = max(coverage_values);
    
    fprintf('Best Coherence: Run %d (%.4f)\n', best_coherence_idx, results.metrics{best_coherence_idx}.max_coherence);
    fprintf('Best Coverage: Run %d (%.2f%%)\n', best_coverage_idx, results.metrics{best_coverage_idx}.column_coverage * 100);
end 