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

% Analyze H-matrix coherence
fprintf('\n--- Analyzing H-Matrix Coherence ---\n');
[max_coherence, ~] = analyze_coherence(H);
fprintf('  Max Coherence: %.4f\n', max_coherence);

%% ===== 6. RECONSTRUCTION =====
fprintf('\n--- Performing Reconstruction ---\n');

% Prepare measurement vector
fprintf('  Preparing measurement vector...\n');
b_vector = echo_data(:);

% Run ADMM reconstruction
fprintf('  Running ADMM reconstruction...\n');
reconstructed_image = run_admm_reconstruction(H, b_vector, params);

%% ===== 7. RESULTS AND VISUALIZATION =====
fprintf('\n--- Results and Visualization ---\n');

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
figure('Position', [100, 100, 800, 600]);
imagesc(imaging_grid.X_mesh * 1000, imaging_grid.Y_mesh * 1000, reconstructed_image);
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

%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix_from_profiles(config, profiles)
    % Initialize Field II
    field_init(-1);
    set_field('fs', config.fs);
    set_field('c', config.c);
    
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
    
    % Channel mapping (matching Arduino CH_PINS array)
    CH_PINS = [26, 24, 28, 29, 30, 31, 32, 33, 34, 36, 37, 39];
    
    % Create element mapping (1-16 for 4x4 array)
    element_mapping = reshape(1:16, num_y_grid, num_x_grid);
    
    % Initialize storage for H-matrix
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
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
        
        % Use quadrant-based selection (as in simulation)
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
        
        % Create apertures
        tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        tx_enabled_matrix(active_tx_indices) = 1;
        TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);
        
        rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
        rx_enabled_matrix(active_rx_indices) = 1;
        RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
        
        % Get profile parameters for this acquisition
        tx_frequencies = profiles.frequencies(r_acq, :);
        delays_us = profiles.delays(r_acq, :);
        apod_weights = profiles.apodizations(r_acq, :);
        tx_phases = profiles.phases(r_acq, :);
        
        % Create excitation signal
        max_len = 0;
        tx_signals = cell(1, config.num_active_tx);
        for k = 1:config.num_active_tx
            f_k = tx_frequencies(k);
            phase_k = tx_phases(k);
            duration = 5 / f_k;
            t_k = 0:1/config.fs:duration;
            signal_k = sin(2*pi*f_k*t_k + phase_k) .* tukeywin(length(t_k), 0.4)';
            tx_signals{k} = signal_k;
            if length(t_k) > max_len, max_len = length(t_k); end
        end
        
        composite_waveform = zeros(1, max_len);
        for k = 1:config.num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
        end
        
        % Set up apertures
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_impulse(RxAperture, ones(1,10));
        xdc_apodization(TxAperture, 0, apod_weights);
        xdc_focus_times(TxAperture, 0, delays_us * 1e-6);
        
        % Calculate impulse response
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        % Clean up apertures
        xdc_free(TxAperture);
        xdc_free(RxAperture);
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    % Assemble H-matrix
    H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, config.fs, 25);
    
    % Clean up Field II
    field_end();
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
    
    % Determine image resolution from H-matrix
    N_pixels = size(H, 2);
    image_size = sqrt(N_pixels);
    if floor(image_size) ~= image_size
        error('H-matrix does not correspond to a square image');
    end
    imageResolution = [image_size, image_size];
    
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