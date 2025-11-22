% TheoryTest2Diag_DiagOnly.m
% Self-contained minimal diagnostics for H and scene/target mapping

addpath(genpath(fullfile(pwd, 'm_files')));

% --- Minimal pipeline to generate H and scene_matrix ---
[output_folder, params] = quantum_init();
params.excitation_amplitude = 1e8; % Diagnostic: very large amplitude
[tx_Aperture, rx_Aperture, imaging_grid] = quantum_pmut_setup(params, output_folder);
H = quantum_h_matrix(tx_Aperture, rx_Aperture, imaging_grid, params, output_folder);
scene_matrix = quantum_scene_creation(imaging_grid, params);

% --- Minimal diagnostics ---
logfile = fullfile('diagnostics_visualized', 'diagnostics_log.txt');
if ~exist('diagnostics_visualized', 'dir')
    mkdir('diagnostics_visualized');
end
logfid = fopen(logfile, 'a');
try
    % H matrix column sums
    col_sums = full(sum(H, 1));
    figure; plot(col_sums); title('Sum of H columns');
    saveas(gcf, fullfile('diagnostics_visualized', 'H_column_sums.png'));

    % v_true_vector and target indices
    v_true_vector = scene_matrix(:);
    nonzero_idx = find(v_true_vector ~= 0);
    col_sums_targets = col_sums(nonzero_idx);
    fprintf('Nonzero indices in v_true_vector: %s\n', mat2str(nonzero_idx'));
    fprintf('Corresponding H column sums: %s\n', mat2str(col_sums_targets));
    fprintf(logfid, 'Nonzero indices in v_true_vector: %s\n', mat2str(nonzero_idx'));
    fprintf(logfid, 'Corresponding H column sums: %s\n', mat2str(col_sums_targets));
    figure; imagesc(reshape(v_true_vector, size(scene_matrix))); colorbar; title('v\_true\_vector reshaped');
    saveas(gcf, fullfile('diagnostics_visualized', 'v_true_vector_reshaped.png'));

    fprintf('Diagnostics and figures saved in %s\n', 'diagnostics_visualized');
    fprintf(logfid, 'Diagnostics and figures saved in %s\n', 'diagnostics_visualized');
catch ME
    fprintf('Diagnostics encountered an error: %s\n', ME.message);
    fprintf(logfid, 'Diagnostics encountered an error: %s\n', ME.message);
end
fclose(logfid);

% --- Local function definitions copied from TheoryTest2Diag ---

function [output_folder, params] = quantum_init()
    % Quantum initialization with JIT compilation and advanced optimizations
    date_str = datestr(now, 'mmddyy');
    base_dir = fullfile(getenv('HOME'), 'Documents', 'MATLAB');
    idx = 1;
    while true
        folder_name = sprintf('%s_%03d', date_str, idx);
        output_folder = fullfile(base_dir, folder_name);
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
            break;
        end
        idx = idx + 1;
    end
    addpath(genpath(fullfile(pwd, 'm_files')));
    field_init(-1);
    feature('jit', 'on');
    feature('accel', 'on');
    maxNumCompThreads('automatic');
    params = struct();
    params.c = 343;
    params.fs = 2e6;
    params.pMUT_width_mm = 20;
    params.pMUT_spacing_mm = 20;
    params.kerf_mm = 0.1;
    params.grid_width_mm = 150;
    params.grid_step_mm = 2;
    params.target_distance_mm = 250;
    params.R_acquisitions = 150;
    params.excitation_amplitude = 1e8;
    params.target_SNR_db = 35;
    params.max_delay_us = 12;
    params.numItersADMM = 50;
    params.rho_admm = 10;
    params.lambda_tv_reg = 0.5;
    params.admm_tol = 1e-7;
    params.admm_max_iter = 50;
    params.pcg_max_iter = 50;
    params.pcg_tol = 1e-10;
    params.use_sparse = true;
    set_field('fs', params.fs);
    set_field('c', params.c);
    quantum_memory_pool = zeros(3000, 3000, 'single');
    clear quantum_memory_pool;
    fprintf('Quantum initialization complete! Output folder: %s\n\n', output_folder);
end

function [tx_Aperture, rx_Aperture, imaging_grid] = quantum_pmut_setup(params, output_folder)
    pMUT_width = params.pMUT_width_mm / 1000;
    pMUT_height = pMUT_width;
    kerf = params.kerf_mm / 1000;
    grid_width = params.grid_width_mm / 1000;
    target_distance = params.target_distance_mm / 1000;
    grid_step = params.grid_step_mm / 1000;
    x_coords_img = -grid_width/2 : grid_step : grid_width/2;
    y_coords_img = -grid_width/2 : grid_step : grid_width/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords_img, y_coords_img);
    N_pixels = numel(X_mesh);
    imaging_grid = struct();
    imaging_grid.x_coords = x_coords_img;
    imaging_grid.y_coords = y_coords_img;
    imaging_grid.z_coord = target_distance;
    imaging_grid.X_mesh = X_mesh;
    imaging_grid.Y_mesh = Y_mesh;
    imaging_grid.N_pixels = N_pixels;
    tx_desired_positions = [
        25e-3, 0, 0;
        -12.5e-3, 21.651e-3, 0;
        -12.5e-3, -21.651e-3, 0
    ];
    rx_pos = [0; 0; 0];
    [tx_Aperture, rx_Aperture] = quantum_aperture_creation(tx_desired_positions, rx_pos, pMUT_width, pMUT_height, kerf, params, output_folder);
end

function [tx_Aperture, rx_Aperture] = quantum_aperture_creation(tx_positions, rx_pos, pMUT_width, pMUT_height, kerf, params, output_folder)
    num_x_grid = 9;
    num_y_grid = 9;
    element_centers = zeros(num_x_grid * num_y_grid, 3);
    center_offset_x = (num_x_grid - 1) / 2 * (pMUT_width + kerf);
    center_offset_y = (num_y_grid - 1) / 2 * (pMUT_height + kerf);
    element_idx = 1;
    for iy = 1:num_y_grid
        y_pos = (iy - 1) * (pMUT_height + kerf) - center_offset_y;
        for ix = 1:num_x_grid
            x_pos = (ix - 1) * (pMUT_width + kerf) - center_offset_x;
            element_centers(element_idx, :) = [x_pos, y_pos, 0];
            element_idx = element_idx + 1;
        end
    end
    tx_active_indices = quantum_element_mapping(tx_positions, element_centers);
    rx_distances = sum((element_centers - repmat(rx_pos', size(element_centers, 1), 1)).^2, 2);
    [~, rx_active_index] = min(rx_distances);
    rx_active_index = rx_active_index(1);
    tx_Aperture = quantum_aperture_matrix(tx_active_indices, num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    rx_Aperture = quantum_aperture_matrix([rx_active_index], num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    quantum_impulse_setup(tx_Aperture, rx_Aperture, params, output_folder);
end

function active_indices = quantum_element_mapping(desired_positions, element_centers)
    active_indices = zeros(size(desired_positions, 1), 1);
    for i = 1:size(desired_positions, 1)
        distances = sum((element_centers - repmat(desired_positions(i, :), size(element_centers, 1), 1)).^2, 2);
        [~, min_idx] = min(distances);
        active_indices(i) = min_idx;
    end
end

function aperture = quantum_aperture_matrix(active_indices, num_x, num_y, width, height, kerf)
    enabled_matrix = zeros(num_y, num_x);
    [row_indices, col_indices] = ind2sub([num_y, num_x], active_indices);
    for i = 1:length(active_indices)
        enabled_matrix(row_indices(i), col_indices(i)) = 1;
    end
    aperture = xdc_2d_array(num_x, num_y, width, height, kerf, kerf, enabled_matrix, 1, 1, [0 0 100e-3]);
end

function quantum_impulse_setup(tx_Aperture, rx_Aperture, params, output_folder)
    fs = params.fs;
    f_min = 45e3;
    f_max = 65e3;
    cycles = 3;
    amplitude = params.excitation_amplitude;
    rng('shuffle');
    tx_frequencies = f_min + (f_max - f_min) * rand(3, 1);
    max_phase_delay = 12e-6;
    tx_phase_delays = max_phase_delay * rand(3, 1);
    tx_durations = cycles ./ tx_frequencies;
    tx_signals = cell(3, 1);
    for i = 1:3
        t_signal = 0 : 1/fs : tx_durations(i);
        signal_base = sin(2 * pi * tx_frequencies(i) * t_signal);
        window = tukeywin(length(signal_base), 0.25)';
        signal_windowed = signal_base .* window;
        tx_signals{i} = signal_windowed * amplitude;
    end
    for i = 1:3
        xdc_impulse(tx_Aperture, tx_signals{i});
        xdc_excitation(tx_Aperture, 1);
    end
    xdc_impulse(rx_Aperture, 1);
    xdc_excitation(rx_Aperture, 1);
end

function H = quantum_h_matrix(tx_Aperture, rx_Aperture, imaging_grid, params, output_folder)
    hydrophone_positions = [imaging_grid.X_mesh(:), imaging_grid.Y_mesh(:), imaging_grid.z_coord * ones(size(imaging_grid.X_mesh(:)))];
    all_hhp_data = cell(params.R_acquisitions, 1);
    all_start_times = zeros(params.R_acquisitions, 1);
    all_K_values = zeros(params.R_acquisitions, 1);
    for r_acq = 1:params.R_acquisitions
        f_min = 45e3;
        f_max = 65e3;
        tx_frequencies = f_min + (f_max - f_min) * rand(3, 1);
        tx_phase_delays = 12e-6 * rand(3, 1);
        cycles = 3;
        tx_durations = cycles ./ tx_frequencies;
        tx_signals = cell(3, 1);
        for i = 1:3
            t_signal = 0 : 1/params.fs : tx_durations(i);
            signal_base = sin(2 * pi * tx_frequencies(i) * t_signal);
            window = tukeywin(length(signal_base), 0.25)';
            tx_signals{i} = signal_base .* window * params.excitation_amplitude;
        end
        for i = 1:3
            xdc_impulse(tx_Aperture, tx_signals{i});
        end
        mean_delay = mean(tx_phase_delays);
        xdc_focus_times(tx_Aperture, 0, [mean_delay, 0, 0]);
        [hhp_r, start_time_r] = calc_hhp(tx_Aperture, rx_Aperture, hydrophone_positions);
        all_hhp_data{r_acq} = hhp_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(hhp_r, 1);
        pause(0.003);
    end
    valid_indices = all_K_values > 0;
    all_end_times = zeros(params.R_acquisitions, 1);
    all_end_times(valid_indices) = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / params.fs;
    min_global_start_time = min(all_start_times);
    max_global_end_time = max(all_end_times);
    if isempty(min_global_start_time) || isempty(max_global_end_time) || min_global_start_time >= max_global_end_time
        min_global_start_time = 0;
        max_K_val = max(all_K_values(all_K_values > 0));
        if isempty(max_K_val) || max_K_val == 0
            max_K_val = 100;
        end
        max_global_end_time = (max_K_val - 1) / params.fs;
        if max_global_end_time < min_global_start_time
            max_global_end_time = min_global_start_time + (100 - 1) / params.fs;
        end
    end
    t_common_axis = min_global_start_time:1/params.fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0
        K_global_per_acq = 1;
        t_common_axis = min_global_start_time;
    end
    total_rows = K_global_per_acq * params.R_acquisitions;
    total_cols = imaging_grid.N_pixels;
    estimated_nnz = total_rows * total_cols * 0.1;
    H_assembled = spalloc(total_rows, total_cols, estimated_nnz);
    current_row_offset = 0;
    for r_acq = 1:params.R_acquisitions
        hhp_current = all_hhp_data{r_acq};
        start_time_current = all_start_times(r_acq);
        K_current = all_K_values(r_acq);
        if K_current == 0 || isempty(hhp_current)
            current_row_offset = current_row_offset + K_global_per_acq;
            continue;
        end
        t_current_acq_axis = start_time_current + (0:(K_current - 1)) / params.fs;
        hhp_aligned_r = quantum_interpolation(t_current_acq_axis, hhp_current, t_common_axis, K_global_per_acq, imaging_grid.N_pixels);
        row_indices = current_row_offset + (1:K_global_per_acq);
        if max(row_indices) <= size(H_assembled, 1)
            H_assembled(row_indices, :) = hhp_aligned_r;
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    H_assembled = sparse(H_assembled);
    H = H_assembled;
end

function hhp_aligned = quantum_interpolation(t_current_acq_axis, hhp_current, t_common_axis, K_global_per_acq, N_pixels)
    hhp_aligned = zeros(K_global_per_acq, N_pixels);
    if length(t_current_acq_axis) > 1 && issorted(t_current_acq_axis)
        for px_col = 1:N_pixels
            if ~isempty(hhp_current) && size(hhp_current, 2) >= px_col
                hhp_aligned(:, px_col) = interp1(t_current_acq_axis, hhp_current(:, px_col), t_common_axis, 'linear', 0);
            end
        end
    elseif isscalar(t_current_acq_axis) && K_global_per_acq >= 1
        [~, idx_match] = min(abs(t_common_axis - t_current_acq_axis));
        if ~isempty(idx_match) && ~isempty(hhp_current)
            hhp_aligned(idx_match, :) = hhp_current(1, :);
        end
    end
end

function scene_matrix = quantum_scene_creation(imaging_grid, params)
    if ~isstruct(imaging_grid) || ~isfield(imaging_grid, 'x_coords') || ~isfield(imaging_grid, 'y_coords')
        warning('imaging_grid is not a struct with x_coords/y_coords. Reconstructing grid.');
        if exist('params','var') && isfield(params, 'grid_width_mm') && isfield(params, 'grid_step_mm')
            grid_width = params.grid_width_mm;
            grid_step = params.grid_step_mm;
            x_coords = -grid_width/2 : grid_step : grid_width/2;
            y_coords = -grid_width/2 : grid_step : grid_width/2;
        else
            x_coords = linspace(-75, 75, 51);
            y_coords = linspace(-75, 75, 51);
        end
    else
        x_coords = imaging_grid.x_coords;
        y_coords = imaging_grid.y_coords;
    end
    scene_matrix = zeros(length(y_coords), length(x_coords));
    grid_positions = [
        -30, 30,  1.0;
        0,   30,  1.0;
        30,  30,  1.0;
        -30, 0,   1.0;
        0,   0,   1.0;
        30,  0,   1.0;
        -30, -30, 1.0;
        0,   -30, 1.0;
        30,  -30, 1.0
    ];
    for i = 1:size(grid_positions, 1)
        x_pos_mm = grid_positions(i, 1);
        y_pos_mm = grid_positions(i, 2);
        amplitude = grid_positions(i, 3);
        [~, ix_scene] = min(abs(x_coords * 1000 - x_pos_mm));
        [~, iy_scene] = min(abs(y_coords * 1000 - y_pos_mm));
        if ix_scene > 0 && ix_scene <= length(x_coords) && ...
           iy_scene > 0 && iy_scene <= length(y_coords)
            scene_matrix(iy_scene, ix_scene) = amplitude;
        end
    end
end 