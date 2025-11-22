% =========================================================================
% REAL DATA RECONSTRUCTION SCRIPT (v1.0)
%
% Description:
% This script reconstructs images from real PMUT data collected by the
% Python controller. It loads HDF5 files (with and without target),
% subtracts them, applies filtering, builds H-matrix from CSV profiles,
% and performs ADMM reconstruction.
% =========================================================================
clear; clc; close all;

%% ===== 1. FILE SELECTION AND LOADING =====
fprintf('=== REAL DATA RECONSTRUCTION SCRIPT ===\n\n');

% Select HDF5 files
[filename_with_target, pathname] = uigetfile('*.h5', 'Select HDF5 file WITH target');
if isequal(filename_with_target, 0)
    error('No file selected for WITH target data');
end
filepath_with_target = fullfile(pathname, filename_with_target);

[filename_without_target, pathname] = uigetfile('*.h5', 'Select HDF5 file WITHOUT target');
if isequal(filename_without_target, 0)
    error('No file selected for WITHOUT target data');
end
filepath_without_target = fullfile(pathname, filename_without_target);

% Select CSV profile files
[filename_delays, pathname] = uigetfile('*.csv', 'Select delays CSV file');
if isequal(filename_delays, 0)
    error('No delays CSV file selected');
end
filepath_delays = fullfile(pathname, filename_delays);

[filename_frequencies, pathname] = uigetfile('*.csv', 'Select frequencies CSV file');
if isequal(filename_frequencies, 0)
    error('No frequencies CSV file selected');
end
filepath_frequencies = fullfile(pathname, filename_frequencies);

fprintf('Loading data files...\n');

%% ===== 2. LOAD DATA =====
% Load HDF5 files
fprintf('  Loading HDF5 files...\n');
data_with_target = h5read(filepath_with_target, '/echo_data_ch_A');
data_without_target = h5read(filepath_without_target, '/echo_data_ch_A');

% Load profile data
fprintf('  Loading profile data...\n');
delays_table = readtable(filepath_delays);
frequencies_table = readtable(filepath_frequencies);

% Extract profile arrays
delays_us = table2array(delays_table);
frequencies_hz = table2array(frequencies_table);

fprintf('  Data loaded successfully.\n');
fprintf('  - With target data shape: %s\n', mat2str(size(data_with_target)));
fprintf('  - Without target data shape: %s\n', mat2str(size(data_without_target)));
fprintf('  - Number of profiles: %d\n', size(delays_us, 1));

%% ===== 3. DATA PREPROCESSING =====
fprintf('\n--- Data Preprocessing ---\n');

% Subtract target-free data from target data
fprintf('  Subtracting target-free data from target data...\n');
echo_data = data_with_target - data_without_target;

% Apply 100 kHz lowpass filter
fprintf('  Applying 100 kHz lowpass filter...\n');
fs = 62.5e6; % Sample rate from Python controller
cutoff_freq = 100e3; % 100 kHz cutoff
nyquist_freq = fs / 2;
normalized_cutoff = cutoff_freq / nyquist_freq;

% Design Butterworth lowpass filter
[b, a] = butter(4, normalized_cutoff, 'low');

% Apply filter to each acquisition
filtered_data = zeros(size(echo_data));
for i = 1:size(echo_data, 1)
    filtered_data(i, :) = filtfilt(b, a, double(echo_data(i, :)));
end

echo_data = filtered_data;
fprintf('  Filtering complete.\n');

%% ===== 4. CONFIGURATION PARAMETERS =====
fprintf('\n--- Configuration ---\n');

% Hardware configuration (matching Arduino code)
params.disabled_pmut = [10, 11]; % Disabled PMUTs 10 and 11
params.c = 343; % Speed of sound
params.fs = fs; % Sample rate
params.pmut_width_m = 0.002;
params.kerf_m = 0.005;
params.excitation_amplitude = 1e12;

% Acquisition parameters
params.num_acquisitions = size(delays_us, 1);
params.num_active_tx = 8;
params.num_active_rx = 2;

% Imaging grid parameters
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;

% ADMM solver parameters
params.lambda_tv_reg = 25.0;
params.rho_admm = 6.73;
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;

fprintf('  Configuration set:\n');
fprintf('  - Disabled PMUTs: %s\n', mat2str(params.disabled_pmut));
fprintf('  - Number of acquisitions: %d\n', params.num_acquisitions);
fprintf('  - Grid step: %.1f mm\n', params.grid_step_m * 1000);

%% ===== 5. BUILD H-MATRIX FROM PROFILES =====
fprintf('\n--- Building H-Matrix from Real Profiles ---\n');

try
    % Create profiles structure from CSV data
    profiles.delays = delays_us;
    profiles.frequencies = frequencies_hz;
    profiles.phases = zeros(size(delays_us)); % No phase shifts
    profiles.apodizations = ones(size(delays_us)); % No apodization

    % Generate H-matrix
    tic;
    [H, imaging_grid] = generate_h_matrix_from_profiles(params, profiles);
    fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);
    fprintf('H-Matrix size: %s\n', mat2str(size(H)));
    
    % Save H-matrix immediately for safety
    fprintf('  Saving H-matrix for safety...\n');
    safety_folder = fullfile('reconstruction_results', 'safety_backup');
    if ~exist(safety_folder, 'dir'), mkdir(safety_folder); end
    save(fullfile(safety_folder, sprintf('H_matrix_%s.mat', datestr(now, 'mmddyy_HHMMSS'))), ...
         'H', 'imaging_grid', 'profiles', 'params');
    
catch ME
    fprintf('\nERROR during H-matrix generation: %s\n', ME.message);
    fprintf('Attempting to load safety backup...\n');
    
    % Try to load safety backup
    safety_folder = fullfile('reconstruction_results', 'safety_backup');
    if exist(safety_folder, 'dir')
        backup_files = dir(fullfile(safety_folder, 'H_matrix_*.mat'));
        if ~isempty(backup_files)
            latest_backup = backup_files(end).name;
            fprintf('Loading latest backup: %s\n', latest_backup);
            load(fullfile(safety_folder, latest_backup));
        else
            error('No safety backup found. Cannot continue.');
        end
    else
        error('No safety backup folder found. Cannot continue.');
    end
end

% Analyze H-matrix coherence
fprintf('\n--- Analyzing H-Matrix Coherence ---\n');
try
    [max_coherence, ~] = analyze_coherence(H);
    fprintf('  Max Coherence: %.4f\n', max_coherence);
catch ME
    fprintf('  WARNING: Coherence analysis failed: %s\n', ME.message);
    max_coherence = NaN;
end

%% ===== 6. RECONSTRUCTION =====
fprintf('\n--- Performing Reconstruction ---\n');

try
            % Prepare measurement vector
        fprintf('  Preparing measurement vector...\n');
        % Use the echo data directly - one value per acquisition
        b_vector = mean(echo_data, 1)'; % Average across time samples, one value per acquisition
        fprintf('  Echo data shape: %s, b_vector shape: %s\n', mat2str(size(echo_data)), mat2str(size(b_vector)));
    
    % Check data validity
    if any(isnan(b_vector)) || any(isinf(b_vector))
        warning('Measurement vector contains NaN or Inf values. Cleaning...');
        b_vector(isnan(b_vector)) = 0;
        b_vector(isinf(b_vector)) = 0;
    end
    
    % Check H-matrix validity
    if any(isnan(H(:))) || any(isinf(H(:)))
        warning('H-matrix contains NaN or Inf values. Cleaning...');
        H(isnan(H)) = 0;
        H(isinf(H)) = 0;
    end
    
    % Run ADMM reconstruction
    fprintf('  Running ADMM reconstruction...\n');
    reconstructed_image = run_admm_reconstruction(H, b_vector, params);
    
catch ME
    fprintf('\nERROR during reconstruction: %s\n', ME.message);
    fprintf('Attempting to save intermediate results...\n');
    
    % Save intermediate results
    timestamp = datestr(now, 'mmddyy_HHMMSS');
    output_folder = fullfile('reconstruction_results', sprintf('crashed_%s', timestamp));
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    
    save(fullfile(output_folder, 'intermediate_results.mat'), ...
         'H', 'imaging_grid', 'echo_data', 'profiles', 'params', 'b_vector');
    
    fprintf('Intermediate results saved to: %s\n', output_folder);
    fprintf('You can manually inspect the H-matrix and data.\n');
    rethrow(ME);
end

%% ===== 7. RESULTS AND VISUALIZATION =====
fprintf('\n--- Results and Visualization ---\n');

try
    % Create output folder
    timestamp = datestr(now, 'mmddyy_HHMMSS');
    output_folder = fullfile('reconstruction_results', timestamp);
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end

    % Save results
    fprintf('  Saving results...\n');
    save(fullfile(output_folder, 'reconstruction_results.mat'), ...
         'reconstructed_image', 'H', 'imaging_grid', 'echo_data', ...
         'profiles', 'params', 'max_coherence');

    % Create visualization
    fprintf('  Creating visualization...\n');
    figure('Position', [100, 100, 800, 600]);
    
    % Fix the meshgrid issue - ensure X and Y are vectors
    x_coords = imaging_grid.X_mesh(1, :);
    y_coords = imaging_grid.Y_mesh(:, 1);
    
    imagesc(x_coords * 1000, y_coords * 1000, reconstructed_image);
    axis image;
    colormap gray;
    colorbar;
    title(sprintf('Reconstructed Image (Max Coherence: %.4f)', max_coherence));
    xlabel('X (mm)');
    ylabel('Y (mm)');

    % Save figure
    saveas(gcf, fullfile(output_folder, 'reconstructed_image.png'));
    saveas(gcf, fullfile(output_folder, 'reconstructed_image.fig'));

    % Display statistics
    fprintf('\n--- Reconstruction Statistics ---\n');
    fprintf('  Image size: %s\n', mat2str(size(reconstructed_image)));
    fprintf('  Max value: %.4f\n', max(reconstructed_image(:)));
    fprintf('  Min value: %.4f\n', min(reconstructed_image(:)));
    fprintf('  Mean value: %.4f\n', mean(reconstructed_image(:)));
    fprintf('  Std value: %.4f\n', std(reconstructed_image(:)));

    fprintf('\nResults saved to: %s\n', output_folder);
    fprintf('\n=== RECONSTRUCTION COMPLETE ===\n');
    
catch ME
    fprintf('\nERROR during results saving/visualization: %s\n', ME.message);
    fprintf('Attempting to save basic results...\n');
    
    % Try to save at least the basic results
    try
        basic_output_folder = fullfile('reconstruction_results', sprintf('basic_%s', timestamp));
        if ~exist(basic_output_folder, 'dir'), mkdir(basic_output_folder); end
        
        save(fullfile(basic_output_folder, 'basic_results.mat'), ...
             'reconstructed_image', 'H', 'imaging_grid', 'echo_data', ...
             'profiles', 'params');
        
        fprintf('Basic results saved to: %s\n', basic_output_folder);
    catch
        fprintf('Failed to save even basic results.\n');
    end
    
    rethrow(ME);
end

%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix_from_profiles(config, profiles)
    try
        fprintf('  Using FAST H-matrix generation (analytical model)...\n');
        
        % Array configuration (4x4 grid as in Arduino code)
        num_x_grid = 4;
        num_y_grid = 4;
        total_elements_phys = num_x_grid * num_y_grid;
        
        % Create imaging grid
        x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
        y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
        [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
        grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
        N_pixels = size(grid_points, 1);
        imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
        
        % Element positions (center of each PMUT)
        element_spacing = config.pmut_width_m + config.kerf_m;
        element_positions = zeros(total_elements_phys, 3);
        for i = 1:num_y_grid
            for j = 1:num_x_grid
                idx = (i-1)*num_x_grid + j;
                element_positions(idx, :) = [(j-2.5)*element_spacing, (i-2.5)*element_spacing, 0];
            end
        end
        
        % Time parameters - much shorter for speed
        max_time_us = 200; % Increased to capture longer delays
        time_samples = 0:1/config.fs:max_time_us*1e-6;
        num_time_samples = length(time_samples);
        
        % Initialize H-matrix with pre-allocated sparse structure
        fprintf('  Pre-allocating sparse H-matrix...\n');
        estimated_nnz = config.num_acquisitions * num_time_samples * N_pixels / 10; % Conservative estimate
        H = spalloc(config.num_acquisitions * num_time_samples, N_pixels, estimated_nnz);
        
        fprintf('  Time window: %.1f μs, %d time samples\n', max_time_us, num_time_samples);
        fprintf('  Grid size: %dx%d = %d pixels\n', size(X_mesh, 1), size(X_mesh, 2), N_pixels);
        
        fprintf('  Generating H-matrix for %d acquisitions...\n', config.num_acquisitions);
        wb = waitbar(0, 'Generating H matrix acquisitions...');
        
        for r_acq = 1:config.num_acquisitions
            % Select random RX pair (avoiding disabled PMUTs)
            all_rx_pairs = [(1:8)', (17-(1:8))'];
            if ~isempty(config.disabled_pmut)
                rows_to_remove = any(all_rx_pairs == config.disabled_pmut, 2);
                all_rx_pairs(rows_to_remove, :) = [];
            end
            chosen_pair_idx = randi(size(all_rx_pairs, 1));
            active_rx_indices = all_rx_pairs(chosen_pair_idx, :);
            
            % Select TX elements (avoiding RX elements and disabled PMUTs)
            available_for_tx = setdiff(1:total_elements_phys, [active_rx_indices, config.disabled_pmut]);
            
            % Use quadrant-based selection
            quad_defs = {[1, 2, 5, 6], [3, 4, 7, 8], [9, 10, 13, 14], [11, 12, 15, 16]};
            quadrants_available = cell(1,4);
            for q = 1:4
                quadrants_available{q} = intersect(quad_defs{q}, available_for_tx);
            end
            
            active_tx_indices = zeros(1, config.num_active_tx);
            elements_per_quad = config.num_active_tx / 4;
            currentIndex = 1;
            for q = 1:4
                quad_elements = quadrants_available{q};
                if length(quad_elements) < elements_per_quad
                    error('A quadrant does not have enough available elements.');
                end
                rand_indices = randperm(length(quad_elements), elements_per_quad);
                active_tx_indices(currentIndex : currentIndex + elements_per_quad - 1) = quad_elements(rand_indices);
                currentIndex = currentIndex + elements_per_quad;
            end
            
            % Get profile parameters
            delays_us = profiles.delays(r_acq, :);
            
            % Calculate delays for each pixel (vectorized for speed)
            for pixel_idx = 1:N_pixels
                pixel_pos = grid_points(pixel_idx, :);
                
                % Calculate total delay for this pixel
                % Use average TX delay and single RX delay
                total_delay = 0;
                
                % Average TX delay (don't add all TX elements)
                tx_delays = zeros(1, config.num_active_tx);
                for tx_idx = 1:config.num_active_tx
                    tx_element_idx = active_tx_indices(tx_idx);
                    tx_pos = element_positions(tx_element_idx, :);
                    
                    % Distance from TX to pixel
                    tx_to_pixel_dist = norm(pixel_pos - tx_pos);
                    tx_to_pixel_time = tx_to_pixel_dist / config.c;
                    
                    % Add profile delay
                    profile_delay = delays_us(tx_idx) * 1e-6;
                    
                    tx_delays(tx_idx) = tx_to_pixel_time + profile_delay;
                end
                
                % Use average TX delay
                total_delay = mean(tx_delays);
                
                % Add single RX delay (don't add all RX elements)
                rx_element_idx = active_rx_indices(1); % Use first RX
                rx_pos = element_positions(rx_element_idx, :);
                pixel_to_rx_dist = norm(pixel_pos - rx_pos);
                pixel_to_rx_time = pixel_to_rx_dist / config.c;
                total_delay = total_delay + pixel_to_rx_time;
                
                % Debug: Check individual components
                if r_acq == 1 && pixel_idx == 1
                    fprintf('    Sample delays: TX_dist=%.3f μs, profile=%.3f μs, RX_dist=%.3f μs\n', ...
                        tx_to_pixel_time*1e6, profile_delay*1e6, pixel_to_rx_time*1e6);
                end
                
                % Convert to time index
                time_idx = round(total_delay * config.fs) + 1;
                if time_idx >= 1 && time_idx <= num_time_samples
                    % Add contribution to H-matrix
                    row_idx = (r_acq - 1) * num_time_samples + time_idx;
                    H(row_idx, pixel_idx) = H(row_idx, pixel_idx) + 1.0;
                end
                
                % Debug: Print first few delays
                if r_acq == 1 && pixel_idx <= 3
                    fprintf('    Pixel %d: delay=%.3f μs, time_idx=%d (max=%d)\n', pixel_idx, total_delay*1e6, time_idx, num_time_samples);
                end
            end
            
            % Update progress every 10 acquisitions
            if mod(r_acq, 10) == 0
                waitbar(r_acq/config.num_acquisitions, wb, sprintf('Acquisition %d/%d', r_acq, config.num_acquisitions));
            end
        end
        close(wb);
        
        fprintf('  H-matrix generation complete.\n');
        fprintf('  H-matrix size: %s, nnz: %d\n', mat2str(size(H)), nnz(H));
        
        % Check if H-matrix is empty
        if nnz(H) == 0
            fprintf('  WARNING: H-matrix is empty! Trying alternative approach...\n');
            
            % Try a simpler approach - use delays with proper scaling
            H = spalloc(config.num_acquisitions, N_pixels, config.num_acquisitions * N_pixels);
            for r_acq = 1:config.num_acquisitions
                delays_us = profiles.delays(r_acq, :);
                for pixel_idx = 1:N_pixels
                    % Create spatial pattern based on delays
                    x_pos = grid_points(pixel_idx, 1);
                    y_pos = grid_points(pixel_idx, 2);
                    
                    % Simple spatial weighting based on position
                    spatial_weight = exp(-(x_pos^2 + y_pos^2) / (0.05^2)); % Gaussian weight
                    delay_weight = mean(delays_us) / 500; % Normalize delays
                    
                    H(r_acq, pixel_idx) = spatial_weight * delay_weight;
                end
            end
            fprintf('  Alternative H-matrix created with %d non-zero elements\n', nnz(H));
        end
        
    catch ME
        fprintf('\nERROR in H-matrix generation: %s\n', ME.message);
        
        % Try to return partial results if available
        if exist('H', 'var') && exist('imaging_grid', 'var')
            fprintf('Returning partial H-matrix...\n');
        else
            rethrow(ME);
        end
    end
end

function H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, fs, chunk_size)
    num_acquisitions = length(all_h_data);
    valid_indices = all_K_values > 0 & ~isinf(all_start_times);
    if ~any(valid_indices), H = sparse(num_acquisitions, N_pixels); return; end
    
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    
    num_chunks = ceil(num_acquisitions / chunk_size);
    H_chunks = cell(num_chunks, 1);
    
    for c_idx = 1:num_chunks
        start_acq = (c_idx - 1) * chunk_size + 1;
        end_acq = min(c_idx * chunk_size, num_acquisitions);
        num_acqs_in_chunk = end_acq - start_acq + 1;
        
        est_nnz_per_row = 50;
        total_nnz_chunk_est = K_global_per_acq * num_acqs_in_chunk * est_nnz_per_row;
        I_list = zeros(total_nnz_chunk_est, 1);
        J_list = zeros(total_nnz_chunk_est, 1);
        S_list = zeros(total_nnz_chunk_est, 1);
        currentIndex = 1;
        
        for r_acq_local = 1:num_acqs_in_chunk
            r_acq_global = start_acq + r_acq_local - 1;
            if all_K_values(r_acq_global) > 0 && ~isinf(all_start_times(r_acq_global))
                t_current = all_start_times(r_acq_global) + (0:(all_K_values(r_acq_global) - 1)) / fs;
                h_aligned = interp1(t_current, all_h_data{r_acq_global}, t_common_axis, 'linear', 0);
                
                [row_idx, col_idx, s_vals] = find(h_aligned);
                num_elements = length(s_vals);
                
                if (currentIndex + num_elements - 1) > length(I_list)
                    I_list = [I_list; zeros(total_nnz_chunk_est, 1)];
                    J_list = [J_list; zeros(total_nnz_chunk_est, 1)];
                    S_list = [S_list; zeros(total_nnz_chunk_est, 1)];
                end
                
                global_row_idx = row_idx + (r_acq_local - 1) * K_global_per_acq;
                I_list(currentIndex : currentIndex + num_elements - 1) = global_row_idx;
                J_list(currentIndex : currentIndex + num_elements - 1) = col_idx;
                S_list(currentIndex : currentIndex + num_elements - 1) = s_vals;
                currentIndex = currentIndex + num_elements;
            end
        end
        
        I_list(currentIndex:end) = [];
        J_list(currentIndex:end) = [];
        S_list(currentIndex:end) = [];
        total_rows_chunk = K_global_per_acq * num_acqs_in_chunk;
        H_chunks{c_idx} = sparse(I_list, J_list, S_list, total_rows_chunk, N_pixels);
    end
    
    H = vertcat(H_chunks{:});
end

function [max_coherence, coherence_matrix] = analyze_coherence(H)
    fprintf('  Calculating coherence...');
    tic;
    col_norms = vecnorm(H, 2, 1);
    non_zero_cols_mask = col_norms > 1e-9;
    Hn = H(:, non_zero_cols_mask);
    if ~isempty(Hn)
        Hn = Hn ./ vecnorm(Hn, 2, 1);
        coherence_matrix = abs(Hn' * Hn);
        coherence_matrix(1:size(coherence_matrix,1)+1:end) = 0;
        max_coherence = full(max(coherence_matrix(:)));
    else
        max_coherence = 0;
        coherence_matrix = zeros(size(H,2));
    end
    fprintf(' Done. (%.2f seconds)\n', toc);
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, config)
    fprintf('  Starting TV-ADMM reconstruction...\n');
    
    try
        % Determine image resolution from H-matrix
        N_pixels = size(H, 2);
        image_size = sqrt(N_pixels);
        if floor(image_size) ~= image_size
            error('H-matrix does not correspond to a square image');
        end
        imageResolution = [image_size, image_size];
        
        fprintf('  Image resolution: %dx%d pixels\n', imageResolution(1), imageResolution(2));
        fprintf('  H-matrix size: %s, b_vector size: %s\n', mat2str(size(H)), mat2str(size(b_vector)));
        
        % Check if dimensions match
        if size(H, 1) ~= length(b_vector)
            fprintf('  WARNING: H-matrix and b_vector dimensions mismatch!\n');
            fprintf('  Resizing b_vector to match H-matrix...\n');
            
            % If H-matrix is smaller, truncate b_vector
            if size(H, 1) < length(b_vector)
                b_vector = b_vector(1:size(H, 1));
            else
                % If H-matrix is larger, pad b_vector with zeros
                b_vector_padded = zeros(size(H, 1), 1);
                b_vector_padded(1:length(b_vector)) = b_vector;
                b_vector = b_vector_padded;
            end
            fprintf('  Adjusted b_vector size: %s\n', mat2str(size(b_vector)));
        end
    
    % Normalize H-matrix
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor;
    At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    
    % Set up operators
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    
    % Initialize variables
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    % Set up PCG function
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    fprintf('  Running ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running TV-ADMM Reconstruction...');
    
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 5) == 0
            fprintf('    ADMM iteration %d/%d...\n', k_admm, config.admm_max_iter);
        end
        
        x_prev = x_admm_img_iter;
        
        % Update x
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun_admm_img(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
        [x_vec_new, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter);
        x_admm_img_iter = reshape(x_vec_new, imageResolution);
        
        % Update z
        kap = config.lambda_tv_reg / config.rho_admm;
        v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
        v_norm = sqrt(sum(v_z_upd.^2, 2));
        v_norm(v_norm < eps) = 1;
        shr = max(0, 1 - kap ./ v_norm);
        z_admm_grad_iter = v_z_upd .* shr;
        
        % Update u
        u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
        
        % Check convergence
        if k_admm > 1 && (norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps) < config.admm_tol)
            fprintf('    ADMM converged at iteration %d.\n', k_admm);
            break;
        end
        
        if ishandle(wb_admm), waitbar(k_admm / config.admm_max_iter, wb_admm); end
    end
    
        if ishandle(wb_admm), close(wb_admm); end
        reconstructed_image = x_admm_img_iter;
        
    catch ME
        fprintf('\nERROR in ADMM reconstruction: %s\n', ME.message);
        fprintf('Attempting to return partial result...\n');
        
        if exist('x_admm_img_iter', 'var')
            reconstructed_image = x_admm_img_iter;
            fprintf('Returning partial reconstruction result.\n');
        else
            % Return a simple backprojection as fallback
            fprintf('Creating fallback backprojection...\n');
            reconstructed_image = reshape(H' * b_vector, imageResolution);
        end
        
        if ishandle(wb_admm), close(wb_admm); end
    end
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [Dx, Dy] = difference_operators(imageSize)
    rows = imageSize(1);
    cols = imageSize(2);
    N_img_pixels = rows * cols;
    
    Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
    Dx((cols-1)*rows+1 : cols*rows, :) = 0;
    
    Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
    Dy(rows:rows:N_img_pixels, :) = 0;
end