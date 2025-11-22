% =========================================================================
% Advanced Comparative H-Matrix Diagnostic Script v3.3
%
% Description:
% This version tests the hypothesis that random transmit apodization will
% create a more incoherent H matrix. It builds and compares three matrices:
%   1. H_hybrid: Driven by the user's REAL experimental data.
%   2. H_simulated: Driven by new random delays and frequencies.
%   3. H_apodized: Driven by new random delays AND random amplitudes.
% =========================================================================

clear; clc; close all;

%% --- 1. User-Defined Parameters & File Selection ---

% --- Data Selection (for H_hybrid) ---
ACQ_START_INDEX = 101;
ACQ_COUNT = 100;

% --- Simulation Parameters ---
c = 343;
fs = 1e6;
pmut_width_m = 0.020;
pmut_spacing_m = 0.020;
kerf_m = 0.0001;
grid_width_m = 0.150;
grid_depth_start_m = 0.250;
grid_depth_end_m = 0.350;
grid_step_m = 0.004;
excitation_amplitude = 1e8;
max_delay_us = 12;

%% --- 2. Setup Output Folder and Load Experimental Data ---
disp('--- Initializing ---');
date_str = datestr(now, 'mmddyy');
base_output_dir = fullfile(pwd, 'diagnostics_output');
if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
counter = 1;
while true
    folder_name = sprintf('comparative_diagnostics_%s_%03d', date_str, counter);
    output_folder = fullfile(base_output_dir, folder_name);
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
        fprintf('Created unique output folder: %s\n', output_folder);
        break;
    end
    counter = counter + 1;
end

[h5_file, h5_path] = uigetfile('*.h5', 'Select the HDF5 Data File');
if isequal(h5_file, 0), disp('User selected Cancel'); return; end
[delays_file, delays_path] = uigetfile('*.csv', 'Select the Delays CSV File', h5_path);
if isequal(delays_file, 0), disp('User selected Cancel'); return; end
[freqs_file, freqs_path] = uigetfile('*.csv', 'Select the Frequencies CSV File', h5_path);
if isequal(freqs_file, 0), disp('User selected Cancel'); return; end

full_delays_us = readmatrix(fullfile(delays_path, delays_file));
full_frequencies_hz = readmatrix(fullfile(freqs_path, freqs_file));

if ACQ_START_INDEX + ACQ_COUNT - 1 > size(full_delays_us, 1)
    error('Selected acquisition range exceeds the number of acquisitions in the file.');
end
acq_indices = ACQ_START_INDEX : (ACQ_START_INDEX + ACQ_COUNT - 1);
exp_delays_us = full_delays_us(acq_indices, :);
exp_frequencies_hz = full_frequencies_hz(acq_indices, :);
num_acquisitions = size(exp_delays_us, 1);
fprintf('Selected %d acquisitions for Hybrid model.\n', num_acquisitions);

%% --- 3. Build All H Matrices ---
disp('--- Building H_hybrid Matrix (experimental parameters) ---');
H_hybrid = generate_h_matrix(num_acquisitions, exp_delays_us, exp_frequencies_hz, c, fs, ...
    pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, ...
    excitation_amplitude, max_delay_us, false); % Apodization OFF

disp('--- Building H_simulated Matrix (random delays & freqs) ---');
rng('default');
sim_delays_us = (max_delay_us * rand(num_acquisitions, 3));
sim_frequencies_hz = (45e3 + (65e3 - 45e3) * rand(num_acquisitions, 3));
H_simulated = generate_h_matrix(num_acquisitions, sim_delays_us, sim_frequencies_hz, c, fs, ...
    pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, ...
    excitation_amplitude, max_delay_us, false); % Apodization OFF

disp('--- Building H_apodized Matrix (random delays & amplitudes) ---');
rng('default'); % Use same random delays for a fair comparison
apodized_delays_us = sim_delays_us;
% For simplicity, use a fixed center frequency when testing apodization
apodized_frequencies_hz = repmat(55000, num_acquisitions, 3);
H_apodized = generate_h_matrix(num_acquisitions, apodized_delays_us, apodized_frequencies_hz, c, fs, ...
    pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, ...
    excitation_amplitude, max_delay_us, true); % Apodization ON

%% --- 4. Run and Save Comparative Diagnostics ---
disp('--- Running and Saving Comparative Diagnostics ---');

x_coords_img = -grid_width_m/2 : grid_step_m : grid_width_m/2;
z_coords_img = grid_depth_start_m : grid_step_m : grid_depth_end_m;
[X_mesh, ~] = meshgrid(x_coords_img, z_coords_img);
imageResolution = size(X_mesh);

run_diagnostics('Hybrid', H_hybrid, fs, num_acquisitions, imageResolution, x_coords_img, z_coords_img, output_folder);
run_diagnostics('Simulated', H_simulated, fs, num_acquisitions, imageResolution, x_coords_img, z_coords_img, output_folder);
run_diagnostics('Apodized', H_apodized, fs, num_acquisitions, imageResolution, x_coords_img, z_coords_img, output_folder);

disp('--- All diagnostics complete! ---');

%% --- Helper Function to Generate an H Matrix ---
function H = generate_h_matrix(num_acquisitions, delays_us, frequencies_hz, c, fs, pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, excitation_amplitude, max_delay_us, use_apodization)
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);

    tx_pos1 = [25e-3, 0, 0];
    tx_pos2 = [-12.5e-3, 21.651e-3, 0];
    tx_pos3 = [-12.5e-3, -21.651e-3, 0];
    tx_desired_positions = [tx_pos1; tx_pos2; tx_pos3];
    rx_pos = [0, 0, 0];
    num_tx_pmuts = size(tx_desired_positions, 1);

    num_x_grid = 9; num_y_grid = 9;
    physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
    center_offset_x = (num_x_grid-1)/2*(pmut_width_m+kerf_m); center_offset_y = (num_y_grid-1)/2*(pmut_width_m+kerf_m);
    for iy = 1:num_y_grid
        y_pos_el = (iy-1)*(pmut_width_m+kerf_m) - center_offset_y;
        for ix = 1:num_x_grid
            x_pos_el = (ix-1)*(pmut_width_m+kerf_m) - center_offset_x;
            element_no = (iy-1)*num_x_grid + ix;
            physical_element_centers(element_no, :) = [x_pos_el, y_pos_el, 0];
        end
    end

    tx_active_indices = zeros(num_tx_pmuts, 1);
    for i = 1:num_tx_pmuts
        distances = sum((physical_element_centers - tx_desired_positions(i, :)).^2, 2);
        [~, min_idx] = min(distances);
        tx_active_indices(i) = min_idx;
    end
    tx_active_indices = unique(tx_active_indices);
    num_tx_active = length(tx_active_indices);

    rx_distances = sum((physical_element_centers - rx_pos).^2, 2);
    [~, rx_active_index] = min(rx_distances);

    tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    tx_enabled_matrix(tx_active_indices) = 1;
    tx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pmut_width_m, pmut_width_m, kerf_m, kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);

    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(rx_active_index) = 1;
    rx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pmut_width_m, pmut_width_m, kerf_m, kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_excitation(rx_Aperture, 1);
    xdc_impulse(rx_Aperture, ones(1,10));

    x_coords_img = -grid_width_m/2 : grid_step_m : grid_width_m/2;
    z_coords_img = grid_depth_start_m : grid_step_m : grid_depth_end_m;
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    
    all_h_data = cell(num_acquisitions, 1);
    all_start_times = zeros(num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix...');
    for r_acq = 1:num_acquisitions
        cycles = 3;
        tx_signals = cell(num_tx_active, 1);
        
        if size(frequencies_hz, 2) == 1
            current_freqs = repmat(frequencies_hz(r_acq, 1), 1, num_tx_active);
        else
            current_freqs = frequencies_hz(r_acq, 1:num_tx_active);
        end
        
        for i = 1:num_tx_active
            tx_durations = cycles / current_freqs(i);
            t_signal = 0 : 1/fs : tx_durations;
            signal_base = sin(2 * pi * current_freqs(i) * t_signal);
            window = tukeywin(length(signal_base), 0.25)';
            tx_signals{i} = signal_base .* window * excitation_amplitude;
        end
        
        xdc_impulse_gorilla(tx_Aperture, tx_signals);
        
        % *** NEW: Apply random apodization if enabled ***
        if use_apodization
            apodization_weights = rand(1, num_tx_active);
            xdc_apodization(tx_Aperture, 0, apodization_weights);
        else
            xdc_apodization(tx_Aperture, 0, ones(1, num_tx_active)); % Default to uniform weights
        end
        
        xdc_excitation(tx_Aperture, ones(1, num_tx_active));
        
        delays_for_sim = delays_us(r_acq, 1:num_tx_active) * 1e-6;
        xdc_focus_times(tx_Aperture, 0, delays_for_sim);
        
        [h_r, start_time_r] = calc_hhp(tx_Aperture, rx_Aperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        waitbar(r_acq/num_acquisitions, wb);
    end
    close(wb);
    
    % Assemble H matrix
    min_start_time = min(all_start_times);
    max_len = 0;
    for i = 1:num_acquisitions
        if ~isempty(all_h_data{i})
            len = size(all_h_data{i}, 1) + round((all_start_times(i) - min_start_time) * fs);
            if len > max_len, max_len = len; end
        end
    end
    K_per_acq = max_len;
    if K_per_acq == 0, K_per_acq = 1; end
    
    H = zeros(K_per_acq * num_acquisitions, N_pixels);
    for r_acq = 1:num_acquisitions
        h_r = all_h_data{r_acq};
        if isempty(h_r), continue; end
        start_time_r = all_start_times(r_acq);
        start_sample_offset = round((start_time_r - min_start_time) * fs);
        num_samples_r = size(h_r, 1);
        
        H_single_acq = zeros(K_per_acq, N_pixels);
        if start_sample_offset > 0 && start_sample_offset + num_samples_r <= K_per_acq + 1
            H_single_acq(start_sample_offset + (1:num_samples_r), :) = h_r;
        end
        
        start_row = (r_acq-1)*K_per_acq + 1;
        end_row = r_acq*K_per_acq;
        H(start_row:end_row, :) = H_single_acq;
    end
    field_end;
end

%% --- Helper Function to Run Diagnostics ---
function run_diagnostics(name, H, fs, num_acquisitions, imageResolution, x_coords_img, z_coords_img, output_folder)
    fprintf('\n--- Running Diagnostics for: %s H-Matrix ---\n', name);
    
    [M, N] = size(H);
    if num_acquisitions == 0, K_per_acq = 0; else, K_per_acq = M / num_acquisitions; end
    
    figure_visibility = 'off';
    log_filepath = fullfile(output_folder, sprintf('diag_%s_log.txt', lower(name)));
    logfid = fopen(log_filepath, 'w');
    fprintf(logfid, '--- Diagnostics for %s H-Matrix ---\n\n', name);
    
    sparsity = nnz(H) / numel(H) * 100;
    fprintf('  H Matrix Size: %d x %d\n', M, N);
    fprintf('  Sparsity: %.4f %% non-zero elements\n', sparsity);
    fprintf(logfid, 'H Matrix Size: %d x %d\n', M, N);
    fprintf(logfid, 'Sparsity: %.4f %% non-zero elements\n', sparsity);
    
    fig1 = figure('Name', sprintf('Diagnostic 1: %s Value Distribution', name), 'visible', figure_visibility);
    histogram(nonzeros(H), 100);
    title(sprintf('Histogram of H Non-Zero Values (%s)', name));
    xlabel('Amplitude'); ylabel('Count'); grid on;
    saveas(fig1, fullfile(output_folder, sprintf('diag_%s_histogram.png', lower(name))));
    close(fig1);

    col_stds = full(std(H,0,1));
    num_zero_cols = sum(col_stds < 1e-20);
    fprintf('  Number of all-zero (or near-zero) columns: %d\n', num_zero_cols);
    fprintf(logfid, 'Number of all-zero (or near-zero) columns: %d\n', num_zero_cols);
    
    max_sv_cols = 500;
    if N > max_sv_cols
        col_subset_indices = randperm(N, max_sv_cols);
        H_subset = H(:, col_subset_indices);
    else
        H_subset = H;
    end
    
    col_norms = vecnorm(H_subset, 2, 1);
    Hn = H_subset ./ col_norms;
    coherence_matrix = abs(Hn' * Hn);
    coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
    max_coherence = max(coherence_matrix(:));
    fprintf('  Max Mutual Coherence: %.4f\n', max_coherence);
    fprintf(logfid, 'Max Mutual Coherence: %.4f\n', max_coherence);

    fprintf('  Calculating RIP proxy (SVD of random submatrices)...\n');
    fprintf(logfid, 'RIP Proxy (SVD of random submatrices):\n');
    K = 10;
    for i = 1:5
        idx = randperm(N, K);
        s = svd(full(H(:,idx)));
        fprintf('    Submatrix %d (K=%d): min(s)=%.3e, max(s)=%.3e\n', i, K, min(s), max(s));
        fprintf(logfid, '    Submatrix %d (K=%d): min(s)=%.3e, max(s)=%.3e\n', i, K, min(s), max(s));
    end
    
    fclose(logfid);
    fprintf('  Diagnostics for %s complete. Log saved.\n', name);
end

% Placeholder for custom multi-impulse function
function xdc_impulse_gorilla(aperture, signals)
    warning('`xdc_impulse_gorilla` is not a standard Field II function. Using the first signal for the entire aperture impulse response.');
    xdc_impulse(aperture, signals{1});
end
