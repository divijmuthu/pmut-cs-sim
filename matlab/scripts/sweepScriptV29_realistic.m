%% sweepScriptV29_realistic.m - Realistic pMUT H-Matrix Optimization
% Based on actual experimental pMUT data: 57.7 kHz resonance, impulse response
% Key insight: pMUT has fixed mechanical resonance, not arbitrary waveforms

clear; clc; close all;

%% Base Configuration (Based on Experimental pMUT Data)
base_config = struct();
base_config.c = 343;                    % Speed of sound (m/s)
base_config.fs = 1e6;                   % Sampling frequency (Hz)
base_config.pmut_width_m = 0.020;       % pMUT width (m)
base_config.tx_pool_width_m = 0.200;    % Transmitter pool width (m) - BEST from V27
base_config.grid_width_m = 0.150;       % Imaging grid width (m)
base_config.target_distance_m = 0.150;  % Target distance (m)
base_config.grid_depth_range_m = 0.100; % Grid depth range (m)
base_config.grid_step_m = 0.010;        % Grid step (m) - BEST from V27
base_config.num_acquisitions = 20;      % BEST from V28_final
base_config.excitation_amplitude = 1e15; % Excitation amplitude

% REALISTIC pMUT PARAMETERS (from experimental data)
base_config.pmut_resonance_freq = 57700; % 57.7 kHz average resonance
base_config.pmut_bandwidth = 2520;      % 2.52 kHz standard deviation
base_config.impulse_duration_us = 10;   % Short impulse excitation (10 μs)

%% Realistic Parameter Sweep (Based on V28_final success)
p = struct();

% BEST performers from V28_final: 5 transmitters, 500μs delays, uniform apodization
p.num_active_tx_sweep = {5};            % Best from V28_final
p.delay_rand_sweep_us = {500};          % Best from V28_final
p.apodization_sweep = {'uniform'};      % Best from V28_final

% Focus on the MOST promising variations from V28_final
p.grid_step_sweep = {0.010};            % BEST grid density from V27
p.tx_pool_width_sweep = {0.200};        % BEST transmitter spacing from V27

% REALISTIC excitation variations (no arbitrary waveforms)
p.impulse_duration_sweep_us = {5, 10, 15}; % Different impulse durations
p.frequency_offset_sweep_hz = {0, 1000, 2000}; % Small offsets from resonance
p.num_acquisitions_sweep = {15, 20, 25}; % Different acquisition counts

%% Output Configuration
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('sweep_output_realistic_v29', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('=== Realistic pMUT H-Matrix Optimization V29 ===\n');
fprintf('Based on actual experimental pMUT data: 57.7 kHz resonance\n');
fprintf('Using realistic impulse excitation (no arbitrary waveforms)\n');
fprintf('Saving results to: %s\n\n', output_folder);

%% Generate Parameter Combinations
param_grid = allcomb(p.num_active_tx_sweep, p.delay_rand_sweep_us, p.apodization_sweep, ...
                     p.grid_step_sweep, p.tx_pool_width_sweep, p.impulse_duration_sweep_us, ...
                     p.frequency_offset_sweep_hz, p.num_acquisitions_sweep);

num_combinations = size(param_grid, 1);
fprintf('Total parameter combinations: %d (realistic pMUT optimization)\n\n', num_combinations);

%% Results Storage (LIGHTWEIGHT)
results = struct();
results.configs = cell(num_combinations, 1);
results.metrics = cell(num_combinations, 1);

%% Main Sweep Loop
for i = 1:num_combinations
    run_name = sprintf('run%03d_tx%d_del%d_%s_gs%.3f_tpw%.3f_imp%d_foff%d_acq%d', i, ...
        param_grid{i, 1}, param_grid{i, 2}, param_grid{i, 3}, ...
        param_grid{i, 4}, param_grid{i, 5}, param_grid{i, 6}, ...
        param_grid{i, 7}, param_grid{i, 8});
    
    fprintf('\n=== [%d/%d] Running Test: %s ===\n', i, num_combinations, run_name);
    fprintf('Parameters: tx=%d, delay=%dμs, apod=%s, grid=%.3fm, pool=%.3fm, impulse=%dμs, foff=%dHz, acq=%d\n', ...
        param_grid{i, 1}, param_grid{i, 2}, param_grid{i, 3}, ...
        param_grid{i, 4}, param_grid{i, 5}, param_grid{i, 6}, ...
        param_grid{i, 7}, param_grid{i, 8});
    
    % Extract parameters
    config = base_config;
    config.num_active_tx = param_grid{i, 1};
    config.max_delay_rand_us = param_grid{i, 2};
    config.apodization_mode = param_grid{i, 3};
    config.grid_step_m = param_grid{i, 4};
    config.tx_pool_width_m = param_grid{i, 5};
    config.impulse_duration_us = param_grid{i, 6};
    config.frequency_offset_hz = param_grid{i, 7};
    config.num_acquisitions = param_grid{i, 8};
    
    fprintf('Starting H matrix generation (realistic pMUT approach)...\n');
    tic;
    
    % Generate H matrix using realistic pMUT behavior
    [H, coherence_matrix] = generate_h_matrix_realistic_pmut(config);
    
    generation_time = toc;
    fprintf('H matrix generation completed in %.2f seconds\n', generation_time);
    
    % Run diagnostics with mutual coherence plot
    fprintf('Running diagnostics with mutual coherence analysis...\n');
    metrics = run_realistic_diagnostics(H, coherence_matrix, config);
    
    % Store results
    results.configs{i} = config;
    results.metrics{i} = metrics;
    
    % Print summary immediately
    fprintf('=== RESULTS for %s ===\n', run_name);
    fprintf('Max Coherence: %.4f (target: <0.85)\n', metrics.max_coherence);
    fprintf('Mean Coherence: %.4f\n', metrics.mean_coherence);
    fprintf('Sparsity: %.2f%% (expected: 75-85%%)\n', metrics.sparsity);
    fprintf('Non-zero Columns: %d/%d\n', full(metrics.num_nonzero_cols), full(metrics.total_cols));
    fprintf('Condition Number: %.2f (target: <100)\n', metrics.condition_number);
    fprintf('RIP Proxy: %.4f (target: <5)\n', metrics.rip_proxy);
    fprintf('Energy Concentration: %.2f%%\n', metrics.energy_concentration);
    fprintf('Rank: %d/%d\n', metrics.rank, min(size(H)));
    
    % Performance assessment
    if metrics.max_coherence < 0.85
        fprintf('🏆 OUTSTANDING: Coherence < 0.85\n');
    elseif metrics.max_coherence < 0.9
        fprintf('✅ EXCELLENT: Coherence < 0.9\n');
    elseif metrics.max_coherence < 0.95
        fprintf('⚠️ GOOD: Coherence < 0.95\n');
    else
        fprintf('❌ POOR: Coherence >= 0.95\n');
    end
    
    if metrics.condition_number < 100
        fprintf('🏆 OUTSTANDING: Condition number < 100\n');
    elseif metrics.condition_number < 200
        fprintf('✅ EXCELLENT: Condition number < 200\n');
    elseif metrics.condition_number < 1000
        fprintf('⚠️ GOOD: Condition number < 1000\n');
    else
        fprintf('❌ POOR: Condition number >= 1000\n');
    end
    
    fprintf('================================\n\n');
    
    % Generate and save mutual coherence plot
    generate_realistic_coherence_plot(coherence_matrix, run_name, output_folder, metrics);
    
    % Clear H matrix to save memory
    clear H coherence_matrix;
end

%% Generate Realistic Summary Report
fprintf('\n=== GENERATING REALISTIC SUMMARY REPORT ===\n');
generate_realistic_summary_report(results, output_folder);

%% Save Complete Results (LIGHTWEIGHT)
save(fullfile(output_folder, 'complete_results.mat'), 'results', 'base_config', 'p');

fprintf('=== Realistic pMUT Optimization Complete! ===\n');
fprintf('Results saved to: %s\n', output_folder);
fprintf('Total disk usage: ~%d MB\n', round(num_combinations * 0.1));

%% Helper Functions

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

function generate_realistic_summary_report(results, output_folder)
    % Create realistic summary CSV
    fprintf('Creating realistic summary CSV...\n');
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
        row.impulse_duration_us = config.impulse_duration_us;
        row.frequency_offset_hz = config.frequency_offset_hz;
        row.num_acquisitions = config.num_acquisitions;
        
        summary_data = [summary_data; row];
    end
    
    % Convert to table and save
    summary_table = struct2table(summary_data);
    writetable(summary_table, fullfile(output_folder, 'realistic_summary.csv'));
    
    % Print best performers with performance assessment
    fprintf('\n=== BEST PERFORMERS (Realistic pMUT) ===\n');
    coherence_values = zeros(length(results.metrics), 1);
    coverage_values = zeros(length(results.metrics), 1);
    condition_values = zeros(length(results.metrics), 1);
    rip_values = zeros(length(results.metrics), 1);
    
    for i = 1:length(results.metrics)
        coherence_values(i) = results.metrics{i}.max_coherence;
        coverage_values(i) = results.metrics{i}.column_coverage;
        condition_values(i) = results.metrics{i}.condition_number;
        rip_values(i) = results.metrics{i}.rip_proxy;
    end
    
    [~, best_coherence_idx] = min(coherence_values);
    [~, best_coverage_idx] = max(coverage_values);
    [~, best_condition_idx] = min(condition_values);
    [~, best_rip_idx] = min(rip_values);
    
    fprintf('Best Coherence: Run %d (%.4f)\n', best_coherence_idx, results.metrics{best_coherence_idx}.max_coherence);
    fprintf('Best Coverage: Run %d (%.2f%%)\n', best_coverage_idx, full(results.metrics{best_coverage_idx}.column_coverage) * 100);
    fprintf('Best Condition Number: Run %d (%.2f)\n', best_condition_idx, results.metrics{best_condition_idx}.condition_number);
    fprintf('Best RIP Proxy: Run %d (%.4f)\n', best_rip_idx, results.metrics{best_rip_idx}.rip_proxy);
    
    % Performance summary
    fprintf('\n=== PERFORMANCE SUMMARY (Realistic pMUT) ===\n');
    outstanding_coherence = sum(coherence_values < 0.85);
    excellent_coherence = sum(coherence_values < 0.9);
    good_coherence = sum(coherence_values < 0.95);
    outstanding_condition = sum(condition_values < 100);
    excellent_condition = sum(condition_values < 200);
    good_condition = sum(condition_values < 1000);
    
    fprintf('Coherence < 0.85: %d/%d runs (%.1f%%) 🏆\n', outstanding_coherence, length(results.metrics), outstanding_coherence/length(results.metrics)*100);
    fprintf('Coherence < 0.9: %d/%d runs (%.1f%%) ✅\n', excellent_coherence, length(results.metrics), excellent_coherence/length(results.metrics)*100);
    fprintf('Coherence < 0.95: %d/%d runs (%.1f%%) ⚠️\n', good_coherence, length(results.metrics), good_coherence/length(results.metrics)*100);
    fprintf('Condition < 100: %d/%d runs (%.1f%%) 🏆\n', outstanding_condition, length(results.metrics), outstanding_condition/length(results.metrics)*100);
    fprintf('Condition < 200: %d/%d runs (%.1f%%) ✅\n', excellent_condition, length(results.metrics), excellent_condition/length(results.metrics)*100);
    fprintf('Condition < 1000: %d/%d runs (%.1f%%) ⚠️\n', good_condition, length(results.metrics), good_condition/length(results.metrics)*100);
    
    % Print all results
    fprintf('\n=== ALL RESULTS (Realistic pMUT) ===\n');
    for i = 1:length(results.metrics)
        fprintf('Run %d: Coherence=%.4f, Coverage=%.2f%%, Cond=%.2f, RIP=%.4f, Impulse=%dμs, Foff=%dHz, Acq=%d\n', ...
            i, results.metrics{i}.max_coherence, results.metrics{i}.column_coverage*100, ...
            results.metrics{i}.condition_number, results.metrics{i}.rip_proxy, ...
            results.configs{i}.impulse_duration_us, results.configs{i}.frequency_offset_hz, results.configs{i}.num_acquisitions);
    end
end 