% =========================================================================
% Advanced Comparative H-Matrix Diagnostic Script v2.1
%
% Description:
% This version fixes a critical crash when building the H_hybrid matrix.
% The script now correctly handles the single-frequency-per-acquisition
% format of the experimental data while still faithfully replicating the
% multi-frequency simulation for the H_simulated matrix.
% =========================================================================

clear; clc; close all;

%% --- 1. User-Defined Parameters & File Selection ---

% --- Data Selection (for H_hybrid) ---
ACQ_START_INDEX = 101;
ACQ_COUNT = 100;

% --- Simulation Parameters (from "Quantum" script) ---
c = 343;                    % Speed of Sound [m/s] (Air)
fs = 1e6;                   % Sampling Frequency [Hz] (Lowered for speed)

% --- PMUT Array Geometry ---
pmut_width_m = 0.020;
pmut_spacing_m = 0.020;
kerf_m = 0.0001;

% --- Imaging Grid Definition ---
grid_width_m = 0.150;
grid_depth_start_m = 0.250;
grid_depth_end_m = 0.350;
grid_step_m = 0.004;

% --- Excitation Parameters ---
excitation_amplitude = 1e8;
max_delay_us = 12;

%% --- 2. Setup Output Folder and Load Experimental Data ---
disp('--- Initializing ---');
% --- Create a unique output folder for this run ---
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

% --- Load Experimental Data for H_hybrid ---
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

%% --- 3. Build H_hybrid (Driven by Experimental Data) ---
disp('--- Building H_hybrid Matrix (using experimental parameters) ---');
H_hybrid = generate_h_matrix(num_acquisitions, exp_delays_us, exp_frequencies_hz, c, fs, ...
    pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, ...
    excitation_amplitude, max_delay_us);

%% --- 4. Build H_simulated (Driven by New Random Data) ---
disp('--- Building H_simulated Matrix (using "Quantum" script logic) ---');
rng('default'); % for reproducibility
sim_delays_us = (max_delay_us * rand(num_acquisitions, 3));
% The simulation uses a different random frequency for each of the 3 transmitters
sim_frequencies_hz = (45e3 + (65e3 - 45e3) * rand(num_acquisitions, 3));
H_simulated = generate_h_matrix(num_acquisitions, sim_delays_us, sim_frequencies_hz, c, fs, ...
    pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, ...
    excitation_amplitude, max_delay_us);

%% --- 5. Run and Save Comparative Diagnostics ---
disp('--- Running and Saving Comparative Diagnostics ---');

x_coords_img = -grid_width_m/2 : grid_step_m : grid_width_m/2;
z_coords_img = grid_depth_start_m : grid_step_m : grid_depth_end_m;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
imageResolution = size(X_mesh);

run_diagnostics('Hybrid', H_hybrid, fs, num_acquisitions, imageResolution, x_coords_img, z_coords_img, output_folder);
run_diagnostics('Simulated', H_simulated, fs, num_acquisitions, imageResolution, x_coords_img, z_coords_img, output_folder);

disp('--- All diagnostics complete! ---');

%% --- Helper Function to Generate an H Matrix ---
function H = generate_h_matrix(num_acquisitions, delays_us, frequencies_hz, c, fs, pmut_width_m, pmut_spacing_m, kerf_m, grid_width_m, grid_depth_start_m, grid_depth_end_m, grid_step_m, excitation_amplitude, max_delay_us)
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
        
        % *** CRITICAL FIX: Handle both single-column and multi-column frequency data ***
        if size(frequencies_hz, 2) == 1
            % Hybrid case: use the same experimental frequency for all transmitters
            current_freqs = repmat(frequencies_hz(r_acq, 1), 1, num_tx_active);
        else
            % Simulated case: use unique random frequencies for each transmitter
            current_freqs = frequencies_hz(r_acq, 1:num_tx_active);
        end
        
        for i = 1:num_tx_active
            tx_durations = cycles / current_freqs(i);
            t_signal = 0 : 1/fs : tx_durations;
            signal_base = sin(2 * pi * current_freqs(i) * t_signal);
            window = tukeywin(length(signal_base), 0.25)';
            tx_signals{i} = signal_base .* window * excitation_amplitude;
        end
        
        % This custom function is needed to set a different impulse for each element
        xdc_impulse_gorilla(tx_Aperture, tx_signals);
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
    if K_per_acq == 0, K_per_acq = 1; end % Prevent zero length
    
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
    
    figure_visibility = 'off'; % Change to 'on' to see plots as they are generated
    
    % --- Diagnostic 1: Energy Map ---
    fprintf('  Calculating energy map...\n');
    column_energies = vecnorm(H, 2, 1);
    energy_map = reshape(column_energies, imageResolution);
    fig1 = figure('Name', sprintf('Diagnostic 1: %s Energy Map', name), 'visible', figure_visibility);
    imagesc(x_coords_img*1000, z_coords_img*1000, energy_map);
    title(sprintf('System Energy Map (%s)', name));
    xlabel('X (mm)'); ylabel('Z (mm)'); colorbar; colormap('hot'); axis image; set(gca, 'YDir', 'normal');
    saveas(fig1, fullfile(output_folder, sprintf('diag_%s_energy_map.png', lower(name))));
    close(fig1);

    % --- Diagnostic 2: Singular Value Decay ---
    fprintf('  Calculating singular values...\n');
    max_sv_cols = 500;
    if N > max_sv_cols
        col_subset_indices = randperm(N, max_sv_cols);
        H_subset = H(:, col_subset_indices);
    else
        H_subset = H;
    end
    S = svd(full(H_subset), 'econ');
    fig2 = figure('Name', sprintf('Diagnostic 2: %s Singular Values', name), 'visible', figure_visibility);
    semilogy(S, 'LineWidth', 1.5);
    title(sprintf('Singular Value Decay (%s)', name));
    xlabel('Singular Value Index'); ylabel('Magnitude'); grid on; axis tight;
    saveas(fig2, fullfile(output_folder, sprintf('diag_%s_singular_values.png', lower(name))));
    close(fig2);
    
    % --- Diagnostic 3: Gram Matrix ---
    fprintf('  Calculating Gram matrix...\n');
    Gram_matrix = H_subset' * H_subset;
    fig3 = figure('Name', sprintf('Diagnostic 3: %s Gram Matrix', name), 'visible', figure_visibility);
    imagesc(abs(Gram_matrix));
    title(sprintf('Gram Matrix |H''*H| (%s)', name));
    xlabel('Pixel Index (Subset)'); ylabel('Pixel Index (Subset)'); colorbar; colormap('hot'); axis square;
    saveas(fig3, fullfile(output_folder, sprintf('diag_%s_gram_matrix.png', lower(name))));
    close(fig3);
    
    % --- Diagnostic 4: Column Time Series ---
    fprintf('  Plotting column time series...\n');
    center_pixel_idx = floor(N / 2);
    h_center_pixel = H(:, center_pixel_idx);
    t_axis_full = (0:length(h_center_pixel)-1) / fs * 1e6;
    fig4 = figure('Name', sprintf('Diagnostic 4: %s Column Time Series', name), 'visible', figure_visibility);
    plot(t_axis_full, h_center_pixel);
    title(sprintf('Spatial Impulse Response from Center Pixel (%s)', name));
    xlabel('Time (us)'); ylabel('Simulated Amplitude'); grid on; axis tight;
    saveas(fig4, fullfile(output_folder, sprintf('diag_%s_column_timeseries.png', lower(name))));
    close(fig4);
    
    fprintf('  Diagnostics for %s complete.\n', name);
end

% This is a placeholder for the custom function from the simulation script.
% It needs to be created in a separate .m file named 'xdc_impulse_gorilla.m'
% on your MATLAB path.
function xdc_impulse_gorilla(aperture, signals)
    % This function is a placeholder. The actual implementation would need to
    % access Field II's internal structures to set per-element impulses,
    % which is not possible with the standard library.
    % For this script, we will just use the first signal for the whole aperture.
    warning('`xdc_impulse_gorilla` is not a standard Field II function. Using the first signal in the cell array for the entire aperture impulse response.');
    xdc_impulse(aperture, signals{1});
end
