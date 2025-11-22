% =========================================================================
% OPTIMAL PROFILE GENERATOR SCRIPT (v1.1 - No Toolboxes)
%
% Description:
% This script implements a greedy algorithm to select a deterministically
% optimized set of acquisition profiles. It generates a large pool of random
% profiles, simulates the acoustic "fingerprint" of each, and then
% iteratively selects the most diverse (least correlated) subset.
%
% v1.1 Improvements:
% - Replaced the toolbox-dependent `corr` function with the standard,
%   built-in `corrcoef` function to ensure compatibility with all MATLAB
%   installations.
% =========================================================================

clear; clc; close all;

%% ===== 1. CONFIGURATION =====
fprintf('=== Optimal Profile Generator (v1.1) ===\n');

% --- Optimization Parameters ---
params.num_profiles_to_generate = 150; % The final number of profiles we want
params.candidate_pool_size = 10000;     % How many random profiles to choose from

% --- Output Setup ---
output_folder = 'optimized_profiles';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving optimized profiles to: %s/\n', output_folder);

% --- Physical Parameters (must match the reconstruction script) ---
params.c = 343;
params.fs = 7.81e6; % Use the downsampled rate for faster simulation (e.g., 31.25e6 / 4)
params.pmut_width_m = 0.020;
params.kerf_m = 0.001;
params.excitation_amplitude = 1e9;
params.num_active_tx = 3;

% --- Profile Randomization Ranges (must match data collection script) ---
params.f_min_hz = 45e3;
params.f_max_hz = 65e3;
params.max_delay_rand_us = 500;
params.max_phase_deg = 360;
params.min_cycles = 1;
params.max_cycles = 4;

%% ===== 2. GENERATE LARGE POOL OF RANDOM PROFILES =====
fprintf('\n--- Step 1: Generating candidate pool of %d random profiles... ---\n', params.candidate_pool_size);
tic;

pool_freqs = randi([params.f_min_hz, params.f_max_hz], params.candidate_pool_size, 1);
pool_delays = randi([0, floor(params.max_delay_rand_us)], params.candidate_pool_size, params.num_active_tx);
pool_phases = randi([0, params.max_phase_deg], params.candidate_pool_size, params.num_active_tx);
pool_cycles = randi([params.min_cycles, params.max_cycles], params.candidate_pool_size, params.num_active_tx);

fprintf('Profile pool generation complete. Time: %.2f seconds.\n', toc);

%% ===== 3. SIMULATE "FINGERPRINT" FOR EACH PROFILE =====
fprintf('\n--- Step 2: Simulating acoustic "fingerprint" for each profile... ---\n');
tic;

% We only need to simulate the response from a single point to get a
% representative fingerprint for correlation analysis.
fingerprint_point = [0, 0, 0.150]; % A single point in the center of the imaging grid

fingerprints = simulate_fingerprints(pool_freqs, pool_delays, pool_phases, pool_cycles, fingerprint_point, params);

fprintf('Fingerprint simulation complete. Time: %.2f seconds.\n', toc);

%% ===== 4. GREEDY SELECTION OF OPTIMAL PROFILES =====
fprintf('\n--- Step 3: Performing greedy selection for %d optimal profiles... ---\n', params.num_profiles_to_generate);
tic;

% --- Initialization ---
pool_indices = 1:params.candidate_pool_size;
optimal_indices = zeros(params.num_profiles_to_generate, 1);
optimal_fingerprints = zeros(size(fingerprints, 1), params.num_profiles_to_generate);

% Start with the first profile
optimal_indices(1) = pool_indices(1);
optimal_fingerprints(:, 1) = fingerprints(:, 1);
pool_indices(1) = []; % Remove from pool

% --- Iteratively select the best next profile ---
wb_greedy = waitbar(0, 'Selecting optimal profiles...');
for i = 2:params.num_profiles_to_generate
    
    candidates_fingerprints = fingerprints(:, pool_indices);
    current_optimal_set = optimal_fingerprints(:, 1:i-1);
    
    % *** FIX: Use the standard `corrcoef` function ***
    % Combine matrices to calculate cross-correlation
    combined_fingerprints = [candidates_fingerprints, current_optimal_set];
    full_corr_matrix = corrcoef(combined_fingerprints);
    
    % Extract the cross-correlation block we care about
    num_candidates = size(candidates_fingerprints, 2);
    corr_matrix = abs(full_corr_matrix(1:num_candidates, num_candidates+1:end));
    
    % Find the max correlation for each candidate
    max_corrs_per_candidate = max(corr_matrix, [], 2);
    
    % Find the candidate with the minimum max correlation
    [~, min_idx_in_pool] = min(max_corrs_per_candidate);
    
    % This is our best new profile
    best_candidate_global_idx = pool_indices(min_idx_in_pool);
    
    % Add it to the optimal set
    optimal_indices(i) = best_candidate_global_idx;
    optimal_fingerprints(:, i) = fingerprints(:, best_candidate_global_idx);
    
    % Remove it from the pool
    pool_indices(min_idx_in_pool) = [];
    
    waitbar(i / params.num_profiles_to_generate, wb_greedy, sprintf('Selecting profile %d/%d', i, params.num_profiles_to_generate));
end
close(wb_greedy);

fprintf('Greedy selection complete. Time: %.2f seconds.\n', toc);

%% ===== 5. SAVE THE OPTIMIZED PROFILES =====
fprintf('\n--- Step 4: Saving optimized profiles to CSV files... ---\n');

% Select the final profiles from the original pools
final_freqs = pool_freqs(optimal_indices, :);
final_delays = pool_delays(optimal_indices, :);
final_phases = pool_phases(optimal_indices, :);
final_cycles = pool_cycles(optimal_indices, :);

% Save to CSV files
writematrix(final_freqs, fullfile(output_folder, 'frequencies_optimized.csv'));
writematrix(final_delays, fullfile(output_folder, 'delays_optimized.csv'));
writematrix(final_phases, fullfile(output_folder, 'phases_optimized.csv'));
writematrix(final_cycles, fullfile(output_folder, 'cycles_optimized.csv'));

fprintf('Successfully saved 4 profile files to the "%s" folder.\n', output_folder);
fprintf('These files can now be used for data collection.\n');

%% ===== HELPER FUNCTION =====

function fingerprints = simulate_fingerprints(freqs, delays, phases, cycles, point, config)
    num_profiles = size(freqs, 1);
    
    % --- Setup Field II (once) ---
    fs = config.fs; c = config.c; num_active_tx = config.num_active_tx;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    num_x_grid = 9; num_y_grid = 9;
    center_offset_x = (num_x_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    center_offset_y = (num_y_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    element_centers = zeros(num_x_grid*num_y_grid, 3);
    for iy = 1:num_y_grid, for ix = 1:num_x_grid
        element_no = (iy-1)*num_x_grid + ix;
        x_pos = (ix-1)*(config.pmut_width_m + config.kerf_m) - center_offset_x;
        y_pos = (iy-1)*(config.pmut_width_m + config.kerf_m) - center_offset_y;
        element_centers(element_no, :) = [x_pos, y_pos, 0];
    end, end
    tx_desired_positions = [ 25e-3, 0, 0; -12.5e-3, 21.651e-3, 0; -12.5e-3, -21.651e-3, 0 ];
    tx_active_indices = zeros(1, num_active_tx);
    for i = 1:num_active_tx, [~, tx_active_indices(i)] = min(sum((element_centers - tx_desired_positions(i,:)).^2, 2)); end
    tx_enabled_matrix = zeros(num_y_grid, num_x_grid); tx_enabled_matrix(tx_active_indices) = 1;
    TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
    [~, rx_active_index] = min(sum(element_centers.^2, 2));
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid); rx_enabled_matrix(rx_active_index) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));

    % --- Simulate each profile ---
    all_h_data = cell(num_profiles, 1);
    all_start_times = zeros(num_profiles, 1);
    all_K_values = zeros(num_profiles, 1);
    wb = waitbar(0, 'Simulating fingerprints...');
    for i = 1:num_profiles
        tx_frequencies = repmat(freqs(i), 1, num_active_tx);
        delays_us = delays(i, :);
        phases_deg = phases(i, :);
        cycle_counts = cycles(i, :);
        
        tx_signals = cell(1, num_active_tx); max_len = 0;
        for k = 1:num_active_tx
            duration = cycle_counts(k) / tx_frequencies(k); t = 0:1/fs:duration;
            phase_rad = phases_deg(k) * pi / 180;
            apodization_factor = cycle_counts(k) / max(cycles(:));
            signal_base = sin(2 * pi * tx_frequencies(k) * t + phase_rad) * apodization_factor;
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx, sig = tx_signals{k}; composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig; end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_excitation(TxAperture, ones(1, num_active_tx));
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, point);
        all_h_data{i} = h_r;
        all_start_times(i) = start_time_r;
        all_K_values(i) = size(h_r, 1);
        waitbar(i/num_profiles, wb);
    end
    close(wb);

    % --- Align fingerprints to a common time axis ---
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    
    fingerprints = zeros(length(t_common_axis), num_profiles);
    for i = 1:num_profiles
        if all_K_values(i) > 0 && ~isempty(all_h_data{i})
            t_current = all_start_times(i) + (0:(all_K_values(i) - 1)) / fs;
            fingerprints(:, i) = interp1(t_current, all_h_data{i}, t_common_axis, 'linear', 0);
        end
    end
    
    xdc_free(TxAperture); xdc_free(RxAperture); field_end();
end
