% =========================================================================
% UNIFIED PMUT SIMULATION SCRIPT (v1.5 - H-Matrix Diagnostics)
%
% Description:
% This script is a diagnostic tool designed to generate and analyze the
% H-matrix from the successful v1.5 simulation. It isolates the H-matrix
% generation and adds a detailed coherence analysis step, including a plot
% of the full coherence matrix. This provides a "gold standard" benchmark
% to compare against the digital twin simulations.fu
%
% =========================================================================
clear; clc; close all;

%% ===== 1. MAIN SIMULATION CONFIGURATION =====
fprintf('=== v1.5 H-Matrix Diagnostic Tool ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('v1_5_diagnostics_output', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s\n\n', output_folder);

% --- Core Parameters from v1.5 ---
params = struct();
params.c = 343;
params.fs = 1e6;
params.pmut_width_m = 0.020;
params.excitation_amplitude = 1e15;

% --- H-Matrix Generation Parameters from v1.5 ---
params.num_acquisitions = 150;
params.num_active_tx = 5;
params.max_delay_rand_us = 500;
params.apodization_mode = 'random';
params.tx_pool_width_m = 0.300;
params.f_min_hz = 45e3;
params.f_max_hz = 65e3;

% --- Imaging Grid Parameters from v1.5 ---
params.grid_width_m = 0.150;
params.target_distance_m = 0.150;
params.grid_depth_range_m = 0.100;
params.grid_step_m = 0.004;


%% ===== 2. GENERATE H-MATRIX =====
fprintf('\n--- Generating H-Matrix using v1.5 logic ---\n');
tic;
[H, imaging_grid] = generate_h_matrix(params);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);


%% ===== 3. ANALYZE H-MATRIX COHERENCE =====
fprintf('\n--- Analyzing H-Matrix Coherence ---\n');
tic;
[max_coherence, coherence_matrix] = analyze_coherence(H);
fprintf('Coherence analysis complete. Time: %.2f seconds.\n', toc);
fprintf('\n  >> Max Coherence of v1.5 H-Matrix: %.4f\n\n', max_coherence);

% --- Visualize Coherence Matrix ---
figure('Visible', 'on', 'Position', [100, 100, 700, 600]);
imagesc(coherence_matrix);
axis square;
colorbar;
title(sprintf('v1.5 Coherence Matrix\nMax Coherence: %.4f', max_coherence));
xlabel('Column Index');
ylabel('Column Index');
saveas(gcf, fullfile(output_folder, 'v1_5_coherence_matrix.png'));

fprintf('Analysis complete. Coherence plot saved to output folder.\n');


%% ===== HELPER FUNCTIONS =====

function [H, imaging_grid] = generate_h_matrix(config)
    fs = config.fs;
    c = config.c;
    num_active_tx = config.num_active_tx;
    
    field_init(-1);
    set_field('fs', fs);
    set_field('c', c);
    
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
    
    % NOTE: The v1.5 sim used an XZ plane. We keep this for a true comparison.
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    grid_points = [X_mesh(:), zeros(numel(X_mesh), 1), Z_mesh(:)];
    N_pixels = size(grid_points, 1);
    
    imaging_grid = struct('X_mesh', X_mesh, 'Z_mesh', Z_mesh, 'points', grid_points);
    all_h_data = cell(config.num_acquisitions, 1);
    all_start_times = zeros(config.num_acquisitions, 1);
    all_K_values = zeros(config.num_acquisitions, 1);
    
    wb = waitbar(0, 'Generating H matrix acquisitions...');
    for r_acq = 1:config.num_acquisitions
        active_indices = randperm(vgrid_total_elements, num_active_tx);
        tx_frequencies = config.f_min_hz + (config.f_max_hz - config.f_min_hz) * rand(1, num_active_tx);
        
        tx_signals = cell(1, num_active_tx);
        max_len = 0;
        for k = 1:num_active_tx
            duration = 3 / tx_frequencies(k);
            t = 0:1/fs:duration;
            random_phase = 2 * pi * rand();
            signal_base = sin(2 * pi * tx_frequencies(k) * t + random_phase);
            window = tukeywin(length(t), 0.25)';
            tx_signals{k} = signal_base .* window;
            if length(t) > max_len, max_len = length(t); end
        end
        composite_waveform = zeros(1, max_len);
        for k = 1:num_active_tx
            sig = tx_signals{k};
            composite_waveform(1:length(sig)) = composite_waveform(1:length(sig)) + sig;
        end
        
        xdc_impulse(TxAperture, composite_waveform * config.excitation_amplitude);
        
        full_apod_vector = zeros(1, vgrid_total_elements);
        apod_weights = ones(1, num_active_tx);
        if strcmp(config.apodization_mode, 'random'), apod_weights = rand(1, num_active_tx); end
        full_apod_vector(active_indices) = apod_weights;
        
        xdc_apodization(TxAperture, 0, full_apod_vector);
        
        full_delay_vector = zeros(1, vgrid_total_elements);
        delays_us = rand(1, num_active_tx) * config.max_delay_rand_us;
        full_delay_vector(active_indices) = delays_us * 1e-6;
        xdc_focus_times(TxAperture, 0, full_delay_vector);
        
        [h_r, start_time_r] = calc_hhp(TxAperture, RxAperture, grid_points);
        all_h_data{r_acq} = h_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(h_r, 1);
        
        waitbar(r_acq/config.num_acquisitions, wb);
    end
    close(wb);
    
    disp('  Assembling H-matrix using interpolation...');
    valid_indices = all_K_values > 0 & all_start_times > -inf;
    if ~any(valid_indices), H = sparse(config.num_acquisitions, N_pixels); return; end
    all_end_times = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / fs;
    min_global_start_time = min(all_start_times(valid_indices));
    max_global_end_time = max(all_end_times);
    t_common_axis = min_global_start_time:1/fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0, K_global_per_acq = 1; end
    total_rows = K_global_per_acq * config.num_acquisitions;
    H = spalloc(total_rows, N_pixels, round(sum(all_K_values) * N_pixels * 0.1));
    current_row_offset = 0;
    for r_acq = 1:config.num_acquisitions
        if all_K_values(r_acq) > 0 && ~isempty(all_h_data{r_acq})
            t_current = all_start_times(r_acq) + (0:(all_K_values(r_acq) - 1)) / fs;
            h_aligned = interp1(t_current, all_h_data{r_acq}, t_common_axis, 'linear', 0);
            row_indices = current_row_offset + (1:K_global_per_acq);
            if max(row_indices) <= size(H, 1)
                H(row_indices, :) = h_aligned;
            end
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    
    xdc_free(TxAperture);
    xdc_free(RxAperture);
    field_end();
end

function [max_coherence, coherence_matrix] = analyze_coherence(H)
    fprintf('  Calculating coherence (this may be slow)...');
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
