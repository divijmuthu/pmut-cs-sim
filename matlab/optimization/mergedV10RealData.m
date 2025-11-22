% =========================================================================
% UNIFIED PMUT ULTRASOUND SIMULATION SCRIPT (v2.1 - Digital Twin)
%
% Description:
% This version has been adapted to be a "digital twin" of the real-world
% experimental setup. It loads experimental data (frequencies, delays, and
% ADC measurements) and uses them to perform a reconstruction.
%
% v2.1 Fixes:
% - Corrected a critical memory error by using spalloc instead of zeros to
%   pre-allocate the large, sparse H-matrix.
% - Clarified the use of the kerf parameter for the virtual grid.
% =========================================================================
clear; clc; close all;

%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== UNIFIED PMUT RECONSTRUCTION FROM REAL DATA ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('real_data_reconstruction_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- File Paths for Real Data (USER MUST EDIT THESE) ---
% ---> ACTION REQUIRED: Enter the full paths to your experimental data files here.
h5_data_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/data_2025-07-18_15-08-13_subtracted.h5';
delays_csv_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/delays_2025-07-18_15-37-31.csv';
freqs_csv_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/frequencies_2025-07-18_15-37-31.csv';

% --- Core Physical and Simulation Parameters ---
params = struct();
params.c = 343;                    % Speed of sound (m/s)
params.fs = 1e6;                   % Sampling frequency (Hz)
params.pmut_width_m = 0.020;       % pMUT width (m)
params.kerf_m = 0.001;             % Kerf for the virtual grid definition
params.excitation_amplitude = 1e15;% High amplitude for strong signal

% --- H-Matrix Generation Parameters (matched to REAL DATA) ---
params.num_acquisitions = 150;     % Using the latter 150 acquisitions
params.num_active_tx = 3;          % HARDWARE CONSTRAINT: 3 transducers
params.use_phase_shifts = false;   % REAL DATA: No random phase shifts

% --- Imaging Grid Parameters ---
params.grid_width_m = 0.150;
params.target_distance_m = 0.150;
params.grid_depth_range_m = 0.100;
params.grid_step_m = 0.004;

% --- ADMM Reconstruction Parameters ---
params.numItersADMM = 50;
params.rho_admm = 6.73;
params.lambda_tv_reg = 1.2;
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;

%% ===== 2. LOAD REAL EXPERIMENTAL DATA =====
fprintf('\n--- Loading Real Experimental Data ---\n');
tic;
% Load delay and frequency profiles
delays_us = readmatrix(delays_csv_path);
frequencies_hz = readmatrix(freqs_csv_path);
% Load raw ADC data
adc_data = h5read(h5_data_path, '/echo_data_raw_adc');

% --- Data Validation and Slicing ---
if size(delays_us, 1) < 200 || size(frequencies_hz, 1) < 200 || size(adc_data, 2) < 200
    error('Data files must contain at least 200 acquisitions.');
end
start_index = size(delays_us, 1) - params.num_acquisitions + 1;
delays_us = delays_us(start_index:end, :);
frequencies_hz = frequencies_hz(start_index:end, :);
adc_data = adc_data(:, start_index:end);

params.TOTAL_SAMPLES = size(adc_data, 1);
fprintf('Loaded %d acquisitions with %d samples each.\n', params.num_acquisitions, params.TOTAL_SAMPLES);
fprintf('Data loading complete. Time: %.2f seconds.\n', toc);

%% ===== 3. GENERATE H-MATRIX FROM EXPERIMENTAL PROFILES =====
fprintf('\n--- Generating H-Matrix from Experimental Profiles ---\n');
tic;
[H, imaging_grid] = generate_h_matrix(params, delays_us, frequencies_hz);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

%% ===== 4. PREPARE MEASUREMENT VECTOR =====
fprintf('\n--- Preparing Measurement Vector from ADC Data ---\n');
b_vector = double(adc_data(:));

%% ===== 5. RECONSTRUCT IMAGE =====
fprintf('\n--- Reconstructing Image via ADMM ---\n');
tic;
reconstructed_image = run_admm_reconstruction(H, b_vector, zeros(size(imaging_grid.X_mesh)), params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);

%% ===== 6. ANALYZE AND SAVE RESULTS =====
fprintf('\n--- Analyzing and Saving Results ---\n');
analyze_and_plot_results(reconstructed_image, H, imaging_grid, params, output_folder);

fprintf('\n\n=== RECONSTRUCTION COMPLETE ===\n');


%% ===== HELPER FUNCTIONS =====

% =========================================================================
% H-MATRIX GENERATION
% =========================================================================
function [H, imaging_grid] = generate_h_matrix(config, delays_us, frequencies_hz)
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);
    
    num_x_grid = 9; num_y_grid = 9;
    total_elements = num_x_grid * num_y_grid;
    center_offset_x = (num_x_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    center_offset_y = (num_y_grid-1)/2*(config.pmut_width_m + config.kerf_m);
    element_centers = zeros(total_elements, 3);
    for iy = 1:num_y_grid, for ix = 1:num_x_grid
        element_no = (iy-1)*num_x_grid + ix;
        x_pos = (ix-1)*(config.pmut_width_m + config.kerf_m) - center_offset_x;
        y_pos = (iy-1)*(config.pmut_width_m + config.kerf_m) - center_offset_y;
        element_centers(element_no, :) = [x_pos, y_pos, 0];
    end, end
    
    tx_desired_positions = [ 25e-3, 0, 0; -12.5e-3, 21.651e-3, 0; -12.5e-3, -21.651e-3, 0 ];
    tx_active_indices = zeros(1, num_active_tx);
    for i = 1:num_active_tx
        distances = sum((element_centers - tx_desired_positions(i,:)).^2, 2);
        [~, tx_active_indices(i)] = min(distances);
    end
    tx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    tx_enabled_matrix(tx_active_indices) = 1;
    TxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0 0 0]);

    [~, rx_active_index] = min(sum(element_centers.^2, 2));
    rx_enabled_matrix = zeros(num_y_grid, num_x_grid);
    rx_enabled_matrix(rx_active_index) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_x_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0 0 0]);
    xdc_impulse(RxAperture, ones(1,10));

    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Z_mesh', Z_mesh, 'points', grid_points);

    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        current_freq = frequencies_hz(r_acq);
        current_delays_us = delays_us(r_acq, :);
        tx_frequencies = repmat(current_freq, 1, num_active_tx);
        
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k);
            t = 0:1/fs:duration;
            random_phase = 0;
            if config.use_phase_shifts, random_phase = 2 * pi * rand(); end
            signal_base = sin(2 * pi * tx_frequencies(k) * t + random_phase);
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end

        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            composite_waveform = composite_waveform + [tx_signals{k}, zeros(1, max_len - length(tx_signals{k}))];
        end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        xdc_excitation(TxAperture, ones(1, num_active_tx));
        xdc_focus_times(TxAperture, 0, current_delays_us * 1e-6);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    disp('  Assembling H-matrix...');
    % *** FIX IS HERE: Use spalloc to prevent memory error ***
    total_rows = config.TOTAL_SAMPLES * config.num_acquisitions;
    estimated_nnz = round(sum(cellfun(@(x) numel(x), all_h_data)) * 0.1); % A better heuristic
    H = spalloc(total_rows, N_pixels, estimated_nnz);
    
    min_start_time = min(all_start_times(all_start_times > -inf));

    for r_acq = 1:config.num_acquisitions
        h_r = all_h_data{r_acq};
        if isempty(h_r), continue; end
        
        start_sample_offset = round((all_start_times(r_acq) - min_start_time) * fs);
        num_samples_r = size(h_r, 1);
        
        start_row_global = (r_acq-1)*config.TOTAL_SAMPLES + 1 + start_sample_offset;
        end_row_global = start_row_global + num_samples_r - 1;
        
        end_row_in_block = min(end_row_global, r_acq*config.TOTAL_SAMPLES);
        num_samples_to_copy = end_row_in_block - start_row_global + 1;

        if num_samples_to_copy > 0
            H(start_row_global:end_row_in_block, :) = h_r(1:num_samples_to_copy, :);
        end
    end

    xdc_free(TxAperture);
    xdc_free(RxAperture);
    field_end();
end

% =========================================================================
% ADMM RECONSTRUCTION
% =========================================================================
function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config)
    imageResolution = size(scene_matrix);
    
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor;
    At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    fprintf('  Starting ADMM iterations...\n');
    for k_admm = 1:config.admm_max_iter
        x_prev = x_admm_img_iter;
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = admm_iteration(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            config.rho_admm, config.lambda_tv_reg, Hfun_pcg_admm, config);
        
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < config.admm_tol, fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        end
    end
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
   rows = imageSize(1); cols = imageSize(2);
   N_img_pixels = rows * cols;
   Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
   Dx( (cols-1)*rows+1 : cols*rows , :) = 0;
   Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
   Dy( rows:rows:N_img_pixels , :) = 0;
end

function [x_new, z_new, u_new] = admm_iteration(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, Hfun_pcg_admm, config)
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm = max(v_norm, eps);
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    u_new = u_old + opDx_tv(x_new) - z_new;
end

% =========================================================================
% ANALYSIS AND VISUALIZATION
% =========================================================================
function analyze_and_plot_results(reconstructed_image, H, imaging_grid, params, output_folder)
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
    end
    fprintf('  Max Coherence of H-Matrix: %.4f\n', max_coherence);
    
    figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
    imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Z_mesh(:,1)*1000, reconstructed_image);
    colormap(gray); axis image; colorbar;
    title('Reconstruction of Real Data');
    xlabel('X Position (mm)'); ylabel('Z Position (mm)');
    saveas(gcf, fullfile(output_folder, 'reconstruction_real_data.png'));
    close(gcf);
    
    metrics = struct('max_coherence', max_coherence);
    save(fullfile(output_folder, 'final_metrics.mat'), 'metrics');
end
