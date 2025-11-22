% =========================================================================
% EXPERIMENTAL DATA RECONSTRUCTION SCRIPT v4.2 (True Digital Twin)
%
% Description:
% This version aligns with the experimental data format from Pythoncontroller.py
% and uses predetermined profiles instead of generating them. It applies a
% 100 kHz lowpass filter and removes phase shifts/apodization for simplicity.
%
% v4.2 Improvements:
% - Aligned with Pythoncontroller.py data format (150 acquisitions, 8 active TX)
% - Added 100 kHz lowpass filter to data preprocessing
% - Removed phase shifts and apodization (uses predetermined profiles)
% - Simplified H-matrix generation to match experimental setup
% =========================================================================
clear; clc; close all;

%% ===== 1. CONFIGURATION =====
fprintf('=== Experimental Data Reconstruction Script (v4.2) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output_real_data', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s/\n', output_folder);

% --- HARDWARE CONFIGURATION (Must match Pythoncontroller.py) ---
params.disabled_pmuts_grid = [10, 11]; 
params.fixed_rx_indices_grid = [3, 14]; 
params.num_active_tx = 8;
params.num_acquisitions = 100; % Updated to match current experiment

% Channel mapping from Pythoncontroller.py
CH_PINS = [26, 24, 28, 29, 30, 31, 32, 33, 34, 36, 37, 39];
% Map grid positions to actual hardware pins
grid_positions = [1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 15, 16];
params.GRID_TO_PIN_MAP = containers.Map(grid_positions, CH_PINS);

% --- Core Physical and Simulation Parameters ---
params.c = 343;                   % Speed of sound in air (m/s)
params.fs = 1e6;                  % Simulation sample rate (Hz)
params.pmut_width_m = 0.002;      % Width of a single pMUT element (m)
params.kerf_m = 0.005;            % Spacing between pMUT elements (m)
params.excitation_amplitude = 1.0;% Amplitude for synthesized waveform

% --- Imaging Grid Parameters ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;

% --- Pre-processing Parameters ---
params.filter_cutoff_hz = 100000;  % 100 kHz lowpass filter
params.filter_order = 4;

% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;
params.assembly_chunk_size = 25;
params.lambda_tv_reg = 25.0;
params.rho_admm = 6.73;

%% ===== 2. LOAD ALL DATA =====
fprintf('\n--- Step 1: Loading All Experimental Data ---\n');
[h5_file, h5_path] = uigetfile('*.h5', 'Select the TARGET HDF5 data file');
if isequal(h5_file, 0), error('User canceled file selection.'); end
data_filepath = fullfile(h5_path, h5_file);
[bg_h5_file, bg_h5_path] = uigetfile('*.h5', 'Select the corresponding BACKGROUND HDF5 file');
if isequal(bg_h5_file, 0), error('Background file is required for reconstruction.'); end
background_filepath = fullfile(bg_h5_path, bg_h5_file);
fprintf('  -> Target: %s\n', data_filepath);
fprintf('  -> Background: %s\n', background_filepath);
fprintf('  Loading datasets from HDF5 files...\n');
target_ch_a_raw = h5read(data_filepath, '/echo_data_ch_A');
target_ch_b_raw = h5read(data_filepath, '/echo_data_ch_B');
background_ch_a_raw = h5read(background_filepath, '/echo_data_ch_A');
background_ch_b_raw = h5read(background_filepath, '/echo_data_ch_B');
tx_pin_profiles = h5read(data_filepath, '/tx_pin_profiles');
fprintf('  Please select the 2 corresponding CSV profile files...\n');
[delays_file, prof_path] = uigetfile('*.csv', 'Select DELAYS_US CSV file');
profiles.delays = readmatrix(fullfile(prof_path, delays_file));
[freqs_file, ~] = uigetfile('*.csv', 'Select FREQUENCIES_HZ CSV file', prof_path);
profiles.frequencies = readmatrix(fullfile(prof_path, freqs_file));

% Automatically set phases to 0 and apodizations to 1
fprintf('  Automatically setting phases to 0 and apodizations to 1...\n');
profiles.phases = zeros(size(profiles.frequencies));
profiles.apodizations = ones(size(profiles.frequencies));
params.num_acquisitions = size(tx_pin_profiles, 2);
fprintf('  Loaded data for %d acquisitions.\n', params.num_acquisitions);

%% ===== 3. PRE-PROCESS EXPERIMENTAL DATA =====
fprintf('\n--- Step 2: Pre-processing Experimental Data ---\n');
% Convert ADC values to voltage (matching Pythoncontroller.py settings)
VOLTAGE_RANGE_MV = 5000.0; % 5V range as per Pythoncontroller.py
RESOLUTION_BITS = 14;
max_adc_value = 2^(RESOLUTION_BITS - 1) - 1;
target_a_mv = (double(target_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
target_b_mv = (double(target_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_a_mv = (double(background_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_b_mv = (double(background_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;

% Get timebase from HDF5 file (matching Pythoncontroller.py)
try, timebase = h5readatt(data_filepath, '/', 'timebase');
catch, fprintf('  INFO: "timebase" attribute not in HDF5 file. Defaulting to 4.\n'); timebase = 4; end
PICO_CLOCK_FREQ_HZ = 62.5e6;
if timebase < 3, sample_interval_s = (2^timebase) / (PICO_CLOCK_FREQ_HZ * 2);
else, sample_interval_s = (timebase - 2) / PICO_CLOCK_FREQ_HZ; end
original_fs = 1 / sample_interval_s;

% Apply 100 kHz lowpass filter
fprintf('  Applying 100 kHz lowpass filter...\n');
[b, a] = butter(params.filter_order, params.filter_cutoff_hz / (original_fs / 2), 'low');
filtered_target_a = filtfilt(b, a, target_a_mv);
filtered_target_b = filtfilt(b, a, target_b_mv);
filtered_background_a = filtfilt(b, a, background_a_mv);
filtered_background_b = filtfilt(b, a, background_b_mv);

fprintf('  Performing background subtraction...\n');
subtracted_a = filtered_target_a - filtered_background_a;
subtracted_b = filtered_target_b - filtered_background_b;
final_a = subtracted_a - mean(subtracted_a, 2);
final_b = subtracted_b - mean(subtracted_b, 2);

fprintf('  Averaging receiver channels...\n');
processed_data = (final_a + final_b) / 2.0;

fprintf('  Resampling data to %.1f MHz...\n', params.fs / 1e6);
[p, q] = rat(params.fs / original_fs);
processed_data_resampled = resample(processed_data', p, q)';
b_vector = processed_data_resampled(:);

%% ===== 4. GENERATE AND ANALYZE H-MATRIX =====
fprintf('\n--- Step 3: Generating Digital Twin H-Matrix ---\n');
tic;
[H_raw, imaging_grid] = generate_h_matrix_experimental(params, profiles, tx_pin_profiles);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

fprintf('  Aligning dimensions of H and b...\n');
num_h_samples = size(H_raw, 1) / params.num_acquisitions;
num_b_samples_raw = length(b_vector) / params.num_acquisitions;
min_samples = floor(min(num_h_samples, num_b_samples_raw));
H_reshaped = reshape(full(H_raw), [floor(num_h_samples), size(H_raw, 2), params.num_acquisitions]);
H_aligned = H_reshaped(1:min_samples, :, :);
H_final = sparse(reshape(H_aligned, [min_samples * params.num_acquisitions, size(H_raw, 2)]));
b_reshaped = reshape(b_vector, [floor(num_b_samples_raw), params.num_acquisitions]);
b_aligned = b_reshaped(1:min_samples, :);
b_final = b_aligned(:);

%% ===== 5. DIAGNOSTIC: COHERENCE ANALYSIS =====
fprintf('\n--- Step 4: Analyzing Final H-Matrix Coherence ---\n');
[max_coherence, coherence_matrix] = analyze_coherence(H_final);
fprintf('  Max Coherence of the final H-matrix: %.4f\n', max_coherence);
figure('Name', 'Coherence Matrix');
imagesc(coherence_matrix);
axis image; colorbar;
title(sprintf('Coherence Matrix (Max Value: %.4f)', max_coherence));
xlabel('Pixel Index'); ylabel('Pixel Index');
fprintf('  Displayed the coherence matrix plot.\n');
fprintf('  A good matrix should be dark everywhere except the main diagonal.\n');
if max_coherence < 1e-6, fprintf('  WARNING: Coherence is zero or near-zero, check simulation.\n'); end
disp('Press any key in the command window to continue to reconstruction...');
pause;

%% ===== 6. RECONSTRUCT IMAGE =====
fprintf('\n--- Step 5: Reconstructing Image via TV-ADMM ---\n');
tic;
scene_matrix_size = size(imaging_grid.X_mesh);
reconstructed_image = run_admm_reconstruction(H_final, b_final, scene_matrix_size, params);
fprintf('ADMM reconstruction complete. Time: %.2f seconds.\n', toc);

%% ===== 7. VISUALIZE AND SAVE =====
fprintf('\n--- Step 6: Visualizing and Saving Results ---\n');
figure('Name', 'Final Reconstruction', 'Position', [100, 100, 700, 600]);
imagesc(imaging_grid.X_mesh(1,:)*1000, imaging_grid.Y_mesh(:,1)*1000, reconstructed_image);
colormap(gray); colorbar; axis image;
title(sprintf('Final Reconstructed Image (%s)', timestamp));
xlabel('X Position (mm)'); ylabel('Y Position (mm)');
saveas(gcf, fullfile(output_folder, 'final_reconstruction.png'));
save(fullfile(output_folder, 'reconstruction_data.mat'), 'reconstructed_image', 'imaging_grid', 'params');
fprintf('Reconstruction saved to: %s\n', output_folder);

%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix_experimental(config, profiles, tx_profiles_pins)
    % This version creates a true digital twin aligned with Pythoncontroller.py
    % experimental setup - uses predetermined profiles, no phase shifts/apodization
    fs = config.fs; c = config.c;
    field_init(-1); set_field('fs', fs); set_field('c', c);
    
    % --- Define geometry (matching Pythoncontroller.py) ---
    num_x_grid = 4; num_y_grid = 4;
    
    % Create reverse mapping from pins to grid positions
    pin_to_grid_map = containers.Map('KeyType','double','ValueType','double');
    all_grid_pos = config.GRID_TO_PIN_MAP.keys;
    all_pins = config.GRID_TO_PIN_MAP.values;
    for i = 1:length(all_pins)
        pin_to_grid_map(all_pins{i}) = all_grid_pos{i};
    end
    
    x_coords_img = -config.grid_x_width_m/2 : config.grid_step_m : config.grid_x_width_m/2;
    y_coords_img = -config.grid_y_width_m/2 : config.grid_step_m : config.grid_y_width_m/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    grid_points = [X_mesh(:), Y_mesh(:), ones(numel(X_mesh), 1) * config.target_height_m];
    N_pixels = size(grid_points, 1);
    imaging_grid = struct('X_mesh', X_mesh, 'Y_mesh', Y_mesh, 'points', grid_points);
    
    % Fixed receiver setup (matching Pythoncontroller.py)
    rx_enabled_matrix = zeros(num_y_grid, num_y_grid);
    rx_enabled_matrix(config.fixed_rx_indices_grid) = 1;
    RxAperture = xdc_2d_array(num_y_grid, num_y_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, rx_enabled_matrix, 1, 1, [0,0,0]);
    xdc_impulse(RxAperture, ones(1,10));
    
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    NUM_MAIN_CYCLES = 7;
    NUM_CANCEL_CYCLES = 1;
    
    wb = waitbar(0, 'Generating Experimental H-Matrix...');
    for r_acq = 1:config.num_acquisitions
        % Get active TX pins from experimental data
        active_tx_pins = tx_profiles_pins(:, r_acq);
        
        % Map pins to grid indices with error checking
        active_grid_indices = [];
        for pin = active_tx_pins
            if isKey(pin_to_grid_map, pin)
                active_grid_indices = [active_grid_indices; pin_to_grid_map(pin)];
            else
                fprintf('Warning: Pin %d not found in mapping for acquisition %d\n', pin, r_acq);
            end
        end
        
        if isempty(active_grid_indices)
            fprintf('Error: No valid grid indices found for acquisition %d\n', r_acq);
            continue;
        end
        
        tx_enabled_matrix = zeros(num_y_grid, num_y_grid);
        tx_enabled_matrix(active_grid_indices) = 1;
        tx_enabled_matrix(config.disabled_pmuts_grid) = 0;
        TxAperture = xdc_2d_array(num_y_grid, num_y_grid, config.pmut_width_m, config.pmut_width_m, config.kerf_m, config.kerf_m, tx_enabled_matrix, 1, 1, [0,0,0]);
        
        % === STEP 1: Synthesize the ideal ANALOG waveform (like Python) ===
        tx_frequencies = profiles.frequencies(r_acq, :);
        tx_phases = profiles.phases(r_acq, :);
        apod_weights = profiles.apodizations(r_acq, :);
        
        % Main Pulse Synthesis
        max_len_main = 0; signals_main = cell(1, config.num_active_tx);
        for k = 1:config.num_active_tx
            duration = NUM_MAIN_CYCLES / tx_frequencies(k); t_k = 0:1/fs:duration;
            signal_k = apod_weights(k) * sin(2*pi*tx_frequencies(k)*t_k + tx_phases(k)) .* blackman(length(t_k))';
            signals_main{k} = signal_k;
            if length(t_k) > max_len_main, max_len_main = length(t_k); end
        end
        composite_main = zeros(1, max_len_main);
        for k = 1:config.num_active_tx, composite_main(1:length(signals_main{k})) = composite_main(1:length(signals_main{k})) + signals_main{k}; end

        % Cancellation Pulse Synthesis
        max_len_cancel = 0; signals_cancel = cell(1, config.num_active_tx);
        for k = 1:config.num_active_tx
            duration = NUM_CANCEL_CYCLES / tx_frequencies(k); t_k = 0:1/fs:duration;
            signal_k = apod_weights(k) * sin(2*pi*tx_frequencies(k)*t_k + tx_phases(k) + pi) .* blackman(length(t_k))';
            signals_cancel{k} = signal_k;
            if length(t_k) > max_len_cancel, max_len_cancel = length(t_k); end
        end
        composite_cancel = zeros(1, max_len_cancel);
        for k = 1:config.num_active_tx, composite_cancel(1:length(signals_cancel{k})) = composite_cancel(1:length(signals_cancel{k})) + signals_cancel{k}; end
        
        analog_waveform = [composite_main, composite_cancel];

        % === STEP 2: Convert to DIGITAL square wave (like Arduino) ===
        zero_crossings = find(diff(sign(analog_waveform)) ~= 0);
        half_period_samples = diff([0, zero_crossings, length(analog_waveform)]);
        digital_waveform = [];
        level = 1;
        for i = 1:length(half_period_samples)
            if sign(analog_waveform(zero_crossings(1))) == 1, level = 1; else, level = -1; end
            digital_waveform = [digital_waveform; ones(half_period_samples(i), 1) * level];
            level = -level;
        end
        digital_waveform = digital_waveform * config.excitation_amplitude;

        % === STEP 3: Set excitation and simulate ===
        xdc_impulse(TxAperture, digital_waveform'); % Use the new digital wave
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        
        % === STEP 4: Apply the single pre-delay (like Arduino) ===
        pre_delay_s = mean(profiles.delays(r_acq, :)) * 1e-6;
        
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r + pre_delay_s; % Add delay here
        all_K_values(r_acq) = size(h_r, 1);
        
        xdc_free(TxAperture);
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    H = assemble_h_matrix_chunked(all_h_data, all_start_times, all_K_values, N_pixels, config.fs, config.assembly_chunk_size);
    xdc_free(RxAperture);
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
        H_chunk_list = cell(num_acqs_in_chunk, 1);
        for r_acq_local = 1:num_acqs_in_chunk
            r_acq_global = start_acq + r_acq_local - 1;
            if all_K_values(r_acq_global) > 0 && ~isinf(all_start_times(r_acq_global))
                t_current = all_start_times(r_acq_global) + (0:(all_K_values(r_acq_global) - 1)) / fs;
                h_aligned = interp1(t_current, all_h_data{r_acq_global}, t_common_axis, 'linear', 0);
                H_chunk_list{r_acq_local} = sparse(h_aligned);
            else
                H_chunk_list{r_acq_local} = sparse(K_global_per_acq, N_pixels);
            end
        end
        H_chunks{c_idx} = vertcat(H_chunk_list{:});
    end
    H = vertcat(H_chunks{:});
end

function reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix_size, config)
    imageResolution = scene_matrix_size;
    H_norm_factor = max(abs(H(:)));
    if H_norm_factor < eps, H_norm_factor = 1; end
    A_admm = H ./ H_norm_factor;
    At_admm = A_admm';
    b_admm_vec = b_vector(:) / H_norm_factor;
    [~, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution);
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + config.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    fprintf('  Starting TV-ADMM iterations...\n');
    wb_admm = waitbar(0, 'Running TV-ADMM Reconstruction...');
    for k_admm = 1:config.admm_max_iter
        if mod(k_admm, 10) == 0, fprintf('  TV-ADMM iteration %d/%d...\n', k_admm, config.admm_max_iter); end
        x_prev = x_admm_img_iter;
        v_upd = z_admm_grad_iter - u_admm_dual_iter;
        bb_upd = Atfun_admm_img(b_admm_vec) + config.rho_admm * opDtx_tv(v_upd);
        [x_vec_new, ~] = pcg(Hfun_pcg_admm, bb_upd(:), config.pcg_tol, config.pcg_max_iter, [], [], x_prev(:));
        x_admm_img_iter = reshape(x_vec_new, imageResolution);
        kap = config.lambda_tv_reg / config.rho_admm;
        v_z_upd = opDx_tv(x_admm_img_iter) + u_admm_dual_iter;
        v_norm = sqrt(sum(v_z_upd.^2, 2)); v_norm(v_norm < eps) = 1;
        shr = max(0, 1 - kap ./ v_norm);
        z_admm_grad_iter = v_z_upd .* shr;
        u_admm_dual_iter = u_admm_dual_iter + opDx_tv(x_admm_img_iter) - z_admm_grad_iter;
        if k_admm > 1 && (norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps) < config.admm_tol), fprintf('    ADMM converged at iteration %d.\n', k_admm); break; end
        if ishandle(wb_admm), waitbar(k_admm / config.admm_max_iter, wb_admm, sprintf('ADMM Iteration %d / %d', k_admm, config.admm_max_iter)); end
    end
    if ishandle(wb_admm), close(wb_admm); end
    reconstructed_image = x_admm_img_iter;
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = operator_setup(A_admm, At_admm, imageResolution)
    Afun_admm = @(x) A_admm * x(:); Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    [Dx_sparse, Dy_sparse] = difference_operators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [Dx, Dy] = difference_operators(imageSize)
   rows = imageSize(1); cols = imageSize(2); N_img_pixels = rows * cols;
   Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
   Dx( (cols-1)*rows+1 : cols*rows , :) = 0;
   Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
   Dy( rows:rows:N_img_pixels , :) = 0;
end