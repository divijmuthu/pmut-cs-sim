%% sweepScriptV26_fast.m - Fast, Focused H-Matrix Optimization
% Based on V23's fast approach with lots of prints and minimal computation
% Focus on the best parameters from V23: 5 transmitters, 500μs delays, uniform apodization

clear; clc; close all;

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
base_config.num_acquisitions = 20;      % Reduced for speed (like V23)
base_config.excitation_amplitude = 1e15; % Excitation amplitude
base_config.f_min_hz = 45000;           % Min frequency (Hz)
base_config.f_max_hz = 65000;           % Max frequency (Hz)

%% Focused Parameter Sweep (Based on V23 Best Results)
p = struct();

% Focus on best performers: 5 transmitters, 500μs delays, uniform apodization
p.num_active_tx_sweep = {5};            % Best from V23
p.delay_rand_sweep_us = {500};          % Best from V23
p.apodization_sweep = {'uniform'};      % Best from V23

% Only test a few key variations
p.grid_step_sweep = {0.008, 0.010, 0.012};  % Grid density optimization
p.tx_pool_width_sweep = {0.200, 0.250};     % Transmitter spacing (avoiding kerf issues)
p.frequency_bandwidth_sweep = {[45000, 65000], [40000, 70000]}; % Frequency diversity

%% Output Configuration
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('sweep_output_fast_v26', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('=== Fast H-Matrix Optimization V26 ===\n');
fprintf('Based on V23 fast approach with focused parameters\n');
fprintf('Saving results to: %s\n\n', output_folder);

%% Generate Parameter Combinations
param_grid = allcomb(p.num_active_tx_sweep, p.delay_rand_sweep_us, p.apodization_sweep, ...
                     p.grid_step_sweep, p.tx_pool_width_sweep, p.frequency_bandwidth_sweep);

num_combinations = size(param_grid, 1);
fprintf('Total parameter combinations: %d (much faster than 81!)\n\n', num_combinations);

%% Results Storage (LIGHTWEIGHT)
results = struct();
results.configs = cell(num_combinations, 1);
results.metrics = cell(num_combinations, 1);

%% Main Sweep Loop
for i = 1:num_combinations
    run_name = sprintf('run%03d_tx%d_del%d_%s_gs%.3f_tpw%.3f_fb%d', i, ...
        param_grid{i, 1}, param_grid{i, 2}, param_grid{i, 3}, ...
        param_grid{i, 4}, param_grid{i, 5}, param_grid{i, 6}(2) - param_grid{i, 6}(1));
    
    fprintf('\n=== [%d/%d] Running Test: %s ===\n', i, num_combinations, run_name);
    fprintf('Parameters: tx=%d, delay=%dμs, apod=%s, grid=%.3fm, pool=%.3fm, freq=%dHz\n', ...
        param_grid{i, 1}, param_grid{i, 2}, param_grid{i, 3}, ...
        param_grid{i, 4}, param_grid{i, 5}, param_grid{i, 6}(2) - param_grid{i, 6}(1));
    
    % Extract parameters
    config = base_config;
    config.num_active_tx = param_grid{i, 1};
    config.max_delay_rand_us = param_grid{i, 2};
    config.apodization_mode = param_grid{i, 3};
    config.grid_step_m = param_grid{i, 4};
    config.tx_pool_width_m = param_grid{i, 5};
    config.f_min_hz = param_grid{i, 6}(1);
    config.f_max_hz = param_grid{i, 6}(2);
    
    fprintf('Starting H matrix generation...\n');
    tic;
    
    % Generate H matrix using V23's fast approach
    H = generate_h_matrix_fast(config);
    
    generation_time = toc;
    fprintf('H matrix generation completed in %.2f seconds\n', generation_time);
    
    % Run fast diagnostics (like V23)
    fprintf('Running diagnostics...\n');
    metrics = run_fast_diagnostics(H, config);
    
    % Store results
    results.configs{i} = config;
    results.metrics{i} = metrics;
    
    % Print summary immediately
    fprintf('=== RESULTS for %s ===\n', run_name);
    fprintf('Max Coherence: %.4f\n', metrics.max_coherence);
    fprintf('Mean Coherence: %.4f\n', metrics.mean_coherence);
    fprintf('Sparsity: %.2f%%\n', metrics.sparsity);
    fprintf('Non-zero Columns: %d/%d\n', metrics.num_nonzero_cols, metrics.total_cols);
    fprintf('Condition Number: %.2f\n', metrics.condition_number);
    fprintf('RIP Proxy: %.4f\n', metrics.rip_proxy);
    fprintf('Energy Concentration: %.2f%%\n', metrics.energy_concentration);
    fprintf('Rank: %d/%d\n', metrics.rank, min(size(H)));
    fprintf('================================\n\n');
    
    % Clear H matrix to save memory
    clear H;
end

%% Generate Summary Report
fprintf('\n=== GENERATING SUMMARY REPORT ===\n');
generate_fast_summary_report(results, output_folder);

%% Save Complete Results (LIGHTWEIGHT)
save(fullfile(output_folder, 'complete_results.mat'), 'results', 'base_config', 'p');

fprintf('=== Fast Optimization Complete! ===\n');
fprintf('Results saved to: %s\n', output_folder);
fprintf('Total disk usage: ~%d MB\n', round(num_combinations * 0.05));

%% Helper Functions

function H = generate_h_matrix_fast(config)
    % Initialize Field II
    fprintf('  Initializing Field II...\n');
    field_init(-1);
    
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
        
        % Randomly choose 3 or 5 transmitters for this acquisition
        if rand < 0.5
            num_active_tx = 3;
        else
            num_active_tx = 5;
        end
        
        fprintf('    Acquisition %d/%d: Using %d transmitters...', r_acq, config.num_acquisitions, num_active_tx);
        
        % Generate diverse excitation signals (like V23)
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        
        % Setup apodization
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random')
            apod_weights = rand(1, num_active_tx);
        end
        
        % Generate composite waveform (like V23)
        excitation_amps = (0.5 + rand(1, num_active_tx)) * config.excitation_amplitude;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k);
            t = 0:1/fs:duration;
            random_phase = 2 * pi * rand();
            signal_base = sin(2 * pi * tx_frequencies(k) * t + random_phase);
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window * excitation_amps(k);
            if length(t) > max_len, max_len = length(t); end
        end
        
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
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
    
    % Assemble H matrix using V23's time-domain approach
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices)
        H = sparse(config.num_acquisitions, N_pixels);
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

function metrics = run_fast_diagnostics(H, config)
    % Fast diagnostics (like V23 but with a few extra metrics)
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
    
    % Coherence analysis (like V23)
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    
    if isempty(Hn) || size(Hn, 2) < 2
        metrics.max_coherence = 0;
        metrics.mean_coherence = 0;
    else
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0; % Remove diagonal
        metrics.max_coherence = full(max(coherence_matrix(:)));
        metrics.mean_coherence = full(mean(coherence_matrix(:)));
        if isempty(metrics.max_coherence), metrics.max_coherence = 0; end
        if isempty(metrics.mean_coherence), metrics.mean_coherence = 0; end
    end
    
    fprintf('    Computing advanced metrics...\n');
    
    % Advanced metrics (simplified)
    if size(Hn, 2) > 1
        metrics.condition_number = cond(full(Hn));
        
        % RIP proxy (simplified)
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

function generate_fast_summary_report(results, output_folder)
    % Create summary CSV
    fprintf('Creating summary CSV...\n');
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
        row.energy_concentration = metrics.energy_concentration;
        row.rank = metrics.rank;
        
        % Configuration parameters
        row.num_active_tx = config.num_active_tx;
        row.max_delay_rand_us = config.max_delay_rand_us;
        row.grid_step_m = config.grid_step_m;
        row.tx_pool_width_m = config.tx_pool_width_m;
        row.f_bandwidth = config.f_max_hz - config.f_min_hz;
        
        summary_data = [summary_data; row];
    end
    
    % Convert to table and save
    summary_table = struct2table(summary_data);
    writetable(summary_table, fullfile(output_folder, 'detailed_summary.csv'));
    
    % Print best performers
    fprintf('\n=== BEST PERFORMERS ===\n');
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
    
    % Print all results
    fprintf('\n=== ALL RESULTS ===\n');
    for i = 1:length(results.metrics)
        fprintf('Run %d: Coherence=%.4f, Coverage=%.2f%%, Cond=%.2f, RIP=%.4f\n', ...
            i, results.metrics{i}.max_coherence, results.metrics{i}.column_coverage*100, ...
            results.metrics{i}.condition_number, results.metrics{i}.rip_proxy);
    end
end 