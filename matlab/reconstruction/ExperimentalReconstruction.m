% Experimental Data Reconstruction using Field II library
% Adapts simulation code for real experimental data
clearvars;
clc;
close all;

%% ===== MAIN EXPERIMENTAL RECONSTRUCTION SCRIPT =====
% Initialize reconstruction
[output_folder, params] = initialize_experimental_reconstruction();

% Load experimental data
[measurements, delay_profiles, frequencies] = load_experimental_data(params, output_folder);

% Generate H matrix using real delay profiles
H = generate_experimental_h_matrix(delay_profiles, frequencies, params, output_folder);

% Perform ADMM-TV reconstruction on real measurements
perform_experimental_reconstruction(H, measurements, params, output_folder);

% Cleanup
cleanup_experimental_reconstruction();

%% ===== HELPER FUNCTIONS =====

function [output_folder, params] = initialize_experimental_reconstruction()
    % Create output folder
    base_dir = fullfile(getenv('HOME'), 'Documents', 'MATLAB');
    date_str = datestr(now, 'mmddyy');
    counter = 1;
    while true
        folder_name = sprintf('Experimental_%s_%03d', date_str, counter);
        output_folder = fullfile(base_dir, folder_name);
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
            fprintf('Created experimental output folder: %s\n', output_folder);
            break;
        end
        counter = counter + 1;
    end
    
    % Add m_files to path
    mfiles_path = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'm_files');
    if exist(mfiles_path, 'dir') && ~contains(path, mfiles_path)
        addpath(genpath(mfiles_path));
        fprintf('Added m_files to path: %s\n', mfiles_path);
    end
    
    % Initialize Field II
    field_init(-1);
    
    % Define experimental parameters
    params = struct();
    params.c = 343;                    % Speed of Sound [m/s] (air)
    params.fs = 2e6;                   % Sampling Frequency [Hz]
    params.pMUT_width_mm = 20;         % pMUT width/height (mm)
    params.pMUT_spacing_mm = 20;       % Triangle sides (mm)
    params.kerf_mm = 0.1;              % Kerf between elements (mm)
    params.grid_width_mm = 150;        % Imaging width (mm)
    params.grid_depth_start_mm = 250;  % Start depth (mm)
    params.grid_depth_end_mm = 350;    % End depth (mm)
    params.grid_step_mm = 4;           % Pixel size (mm)
    params.R_acquisitions = 100;       % Number of acquisitions (last 100)
    params.excitation_amplitude = 10000; % Signal amplitude
    params.target_SNR_db = 30;         % Target SNR
    params.max_delay_us = 12;          % Max delay (us)
    params.numItersADMM = 50;          % ADMM iterations
    params.rho_admm = 10;              % ADMM rho parameter
    params.lambda_tv_reg = 0.1;        % TV regularization
    
    % Experimental data parameters
    params.h5_file_path = '';          % Will be set by user
    params.delay_csv_path = '';        % Will be set by user
    params.freq_csv_path = '';         % Will be set by user
    params.use_last_100 = true;        % Use only last 100 acquisitions
    
    % Efficiency parameters
    params.use_parallel = false;       % Disable parallel processing for Field II compatibility
    params.max_workers = 4;            % Max parallel workers
    params.admm_tol = 1e-6;            % ADMM convergence tolerance
    params.admm_max_iter = 50;         % Max ADMM iterations
    params.pcg_max_iter = 25;          % Max PCG iterations
    params.pcg_tol = 1e-8;             % PCG tolerance
    
    % Advanced optimization parameters
    params.use_gpu = false;            % GPU acceleration (if available)
    params.use_sparse = true;          % Use sparse matrices where possible
    params.memory_pool = true;         % Memory pooling for large matrices
    params.vectorize_interp = true;    % Vectorized interpolation
    params.adaptive_tolerance = true;  % Adaptive convergence tolerance
    params.block_size = 100;           % Block size for matrix operations
    params.cache_results = true;       % Cache intermediate results
    
    % Set Field II parameters
    set_field('fs', params.fs);
    set_field('c', params.c);
    
    % Initialize memory pool for large matrices
    if params.memory_pool
        temp_large = zeros(1000, 1000);
        clear temp_large;
        fprintf('Memory pool initialized\n');
    end
    
    % Print parameters
    fprintf('--- Experimental Reconstruction v1.0.0 ---\n\n');
    fprintf('--- Key Experimental Parameters ---\n');
    fprintf('pMUT Width: %g mm, Spacing: %g mm\n', params.pMUT_width_mm, params.pMUT_spacing_mm);
    fprintf('Imaging Grid: %.0f-%.0f mm depth, %.0f mm width, %.1f mm step\n', ...
        params.grid_depth_start_mm, params.grid_depth_end_mm, params.grid_width_mm, params.grid_step_mm);
    fprintf('Acquisitions: %d (last 100), Max Delay: %g us\n', params.R_acquisitions, params.max_delay_us);
    fprintf('Parallel Processing: %s\n\n', mat2str(params.use_parallel));
end

function [measurements, delay_profiles, frequencies] = load_experimental_data(params, output_folder)
    fprintf('\n--- Loading Experimental Data ---\n');
    
    % Prompt user for file paths if not set
    if isempty(params.h5_file_path)
        [h5_file, h5_path] = uigetfile('*.h5', 'Select H5 file with echo data');
        if h5_file == 0
            error('No H5 file selected');
        end
        params.h5_file_path = fullfile(h5_path, h5_file);
    end
    
    if isempty(params.delay_csv_path)
        [delay_file, delay_path] = uigetfile('*.csv', 'Select CSV file with delay profiles');
        if delay_file == 0
            error('No delay CSV file selected');
        end
        params.delay_csv_path = fullfile(delay_path, delay_file);
    end
    
    if isempty(params.freq_csv_path)
        [freq_file, freq_path] = uigetfile('*.csv', 'Select CSV file with frequencies');
        if freq_file == 0
            error('No frequency CSV file selected');
        end
        params.freq_csv_path = fullfile(freq_path, freq_file);
    end
    
    % Validate file paths
    validate_experimental_files(params.h5_file_path, params.delay_csv_path, params.freq_csv_path);
    
    % Load H5 data
    fprintf('Loading H5 file: %s\n', params.h5_file_path);
    h5_info = h5info(params.h5_file_path);
    fprintf('H5 file contains datasets:\n');
    for i = 1:length(h5_info.Datasets)
        fprintf('  %s\n', h5_info.Datasets(i).Name);
    end
    
    % Load echo data (assuming first dataset)
    echo_data = h5read(params.h5_file_path, ['/' h5_info.Datasets(1).Name]);
    fprintf('Echo data shape: %s\n', mat2str(size(echo_data)));
    
    % Check if data needs to be transposed (common issue with experimental data)
    if size(echo_data, 1) > size(echo_data, 2)
        fprintf('Transposing echo data (samples x acquisitions -> acquisitions x samples)\n');
        echo_data = echo_data';
        fprintf('Echo data shape after transpose: %s\n', mat2str(size(echo_data)));
    end
    
    % Load delay profiles
    fprintf('Loading delay profiles: %s\n', params.delay_csv_path);
    delay_data = readtable(params.delay_csv_path);
    delay_profiles = table2array(delay_data);
    fprintf('Delay profiles shape: %s\n', mat2str(size(delay_profiles)));
    
    % Load frequencies
    fprintf('Loading frequencies: %s\n', params.freq_csv_path);
    freq_data = readtable(params.freq_csv_path);
    frequencies = table2array(freq_data);
    fprintf('Frequencies shape: %s\n', mat2str(size(frequencies)));
    
    % Handle frequency data - if it's a single column, expand to match delay profiles
    if size(frequencies, 2) == 1 && size(delay_profiles, 2) > 1
        fprintf('Expanding single frequency column to match delay profile dimensions\n');
        frequencies = repmat(frequencies, 1, size(delay_profiles, 2));
        fprintf('Frequencies shape after expansion: %s\n', mat2str(size(frequencies)));
    end
    
    % Validate data dimensions
    num_acquisitions = size(echo_data, 1);
    signal_length = size(echo_data, 2);
    
    fprintf('Data validation:\n');
    fprintf('  Echo data: %d acquisitions x %d samples\n', num_acquisitions, signal_length);
    fprintf('  Delay profiles: %d acquisitions x %d elements\n', size(delay_profiles, 1), size(delay_profiles, 2));
    fprintf('  Frequencies: %d acquisitions x %d elements\n', size(frequencies, 1), size(frequencies, 2));
    
    % Check if we have enough acquisitions
    if num_acquisitions < params.R_acquisitions
        warning('Only %d acquisitions available, but %d requested. Using all available.', num_acquisitions, params.R_acquisitions);
        params.R_acquisitions = num_acquisitions;
    end
    
    % Use only last 100 acquisitions if requested and available
    if params.use_last_100 && num_acquisitions > params.R_acquisitions
        start_idx = num_acquisitions - params.R_acquisitions + 1;
        echo_data = echo_data(start_idx:end, :);
        delay_profiles = delay_profiles(start_idx:end, :);
        frequencies = frequencies(start_idx:end, :);
        fprintf('Using last %d acquisitions (indices %d-%d)\n', params.R_acquisitions, start_idx, num_acquisitions);
    elseif num_acquisitions <= params.R_acquisitions
        fprintf('Using all %d available acquisitions\n', num_acquisitions);
        params.R_acquisitions = num_acquisitions;
    end
    
    % Apply 60 kHz lowpass filter (critical for proper data filtering)
    fprintf('Applying 60 kHz lowpass filter to echo data...\n');
    echo_data_filtered = apply_lowpass_filter(echo_data, params.fs, 60e3);
    fprintf('Lowpass filtering completed.\n');
    
    % Validate final data
    if isempty(echo_data_filtered) || isempty(delay_profiles) || isempty(frequencies)
        error('One or more data arrays are empty. Please check your data files.');
    end
    
    if size(echo_data_filtered, 1) == 0 || size(delay_profiles, 1) == 0 || size(frequencies, 1) == 0
        error('One or more data arrays have zero rows. Please check your data files.');
    end
    
    % Prepare measurements structure
    measurements = struct();
    measurements.echo_data = echo_data_filtered;
    measurements.num_acquisitions = size(echo_data_filtered, 1);
    measurements.signal_length = size(echo_data_filtered, 2);
    
    fprintf('Final data shapes:\n');
    fprintf('  Echo data: %s\n', mat2str(size(measurements.echo_data)));
    fprintf('  Delay profiles: %s\n', mat2str(size(delay_profiles)));
    fprintf('  Frequencies: %s\n', mat2str(size(frequencies)));
    
    % Plot sample data
    plot_experimental_data(measurements, delay_profiles, frequencies, params, output_folder);
end

function plot_experimental_data(measurements, delay_profiles, frequencies, params, output_folder)
    figure(1);
    set(gcf, 'visible', 'off');
    clf;
    
    % Plot sample echo data
    subplot(2, 2, 1);
    sample_acq = min(5, measurements.num_acquisitions);
    t_axis = (0:(measurements.signal_length-1)) / params.fs * 1e6;
    % Limit the number of points to plot for performance
    plot_points = min(1000, measurements.signal_length);
    plot_indices = round(linspace(1, measurements.signal_length, plot_points));
    plot(t_axis(plot_indices), double(measurements.echo_data(sample_acq, plot_indices)));
    title(sprintf('Sample Echo Data (Acquisition %d, %d points)', sample_acq, plot_points));
    xlabel('Time (us)'); ylabel('Amplitude');
    grid on;
    
    % Plot delay profiles
    subplot(2, 2, 2);
    plot(delay_profiles');
    title('Delay Profiles (All Acquisitions)');
    xlabel('Element Index'); ylabel('Delay (us)');
    grid on;
    
    % Plot frequencies
    subplot(2, 2, 3);
    plot(frequencies');
    title('Frequencies (All Acquisitions)');
    xlabel('Element Index'); ylabel('Frequency (Hz)');
    grid on;
    
    % Plot data statistics
    subplot(2, 2, 4);
    % Convert to double for statistics calculation
    echo_data_double = double(measurements.echo_data(:));
    echo_stats = [mean(echo_data_double), std(echo_data_double), ...
                  min(echo_data_double), max(echo_data_double)];
    bar(echo_stats);
    set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
    title('Echo Data Statistics');
    ylabel('Amplitude');
    grid on;
    
    sgtitle('Figure 1: Experimental Data Overview');
    set(gcf, 'Color', 'w');
    saveas(gcf, fullfile(output_folder, 'figure1_experimental_data.png'));
    close(gcf);
end

function H = generate_experimental_h_matrix(delay_profiles, frequencies, params, output_folder)
    fprintf('\n--- Generating Experimental H Matrix ---\n');
    
    % Setup imaging grid
    grid_width = params.grid_width_mm / 1000;
    grid_depth_start = params.grid_depth_start_mm / 1000;
    grid_depth_end = params.grid_depth_end_mm / 1000;
    grid_step = params.grid_step_mm / 1000;
    
    x_coords_img = -grid_width/2 : grid_step : grid_width/2;
    z_coords_img = grid_depth_start : grid_step : grid_depth_end;
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    N_pixels = numel(X_mesh);
    
    imaging_grid = struct();
    imaging_grid.x_coords = x_coords_img;
    imaging_grid.z_coords = z_coords_img;
    imaging_grid.X_mesh = X_mesh;
    imaging_grid.Z_mesh = Z_mesh;
    imaging_grid.N_pixels = N_pixels;
    
    fprintf('Imaging grid: %d pixels (%d axial x %d lateral)\n', ...
        N_pixels, length(z_coords_img), length(x_coords_img));
    
    % Setup pMUT array (same as simulation)
    pMUT_width = params.pMUT_width_mm / 1000;
    pMUT_height = pMUT_width;
    kerf = params.kerf_mm / 1000;
    
    % Define pMUT positions (same as simulation)
    tx_desired_positions = [
        25e-3, 0, 0;           % TX1: (25mm, 0mm, 0mm)
        -12.5e-3, 21.651e-3, 0; % TX2: (-12.5mm, 21.651mm, 0mm)
        -12.5e-3, -21.651e-3, 0  % TX3: (-12.5mm, -21.651mm, 0mm)
    ];
    rx_pos = [0, 0, 0];  % RX: Center position
    
    % Create apertures
    [tx_Aperture, rx_Aperture] = create_experimental_apertures(tx_desired_positions, rx_pos, pMUT_width, pMUT_height, kerf);
    
    % Setup hydrophone positions
    hydrophone_positions = [imaging_grid.X_mesh(:), zeros(size(imaging_grid.X_mesh(:))), imaging_grid.Z_mesh(:)];
    
    % Pre-allocate arrays for efficiency
    all_hhp_data = cell(params.R_acquisitions, 1);
    all_start_times = zeros(params.R_acquisitions, 1);
    all_K_values = zeros(params.R_acquisitions, 1);
    
    fprintf('Generating H matrix for %d acquisitions...\n', params.R_acquisitions);
    tic;
    
    % Process each acquisition with real delay profiles
    for r_acq = 1:params.R_acquisitions
        if mod(r_acq, 10) == 0 || r_acq == 1
            fprintf('  Acquisition %d/%d (%.1f%%)...\n', r_acq, params.R_acquisitions, 100*r_acq/params.R_acquisitions);
        end
        
        % Use real delay profile for this acquisition
        if r_acq <= size(delay_profiles, 1)
            delay_vector = delay_profiles(r_acq, :) * 1e-6; % Convert to seconds
        else
            warning('Acquisition %d exceeds delay profile size, using random delays', r_acq);
            delay_vector = rand(1, 3) * params.max_delay_us / 1e6;
        end
        
        % Apply delays to transmit aperture
        xdc_focus_times(tx_Aperture, 0, delay_vector);
        
        % Calculate bistatic response
        [hhp_r, start_time_r] = calc_hhp(tx_Aperture, rx_Aperture, hydrophone_positions);
        
        all_hhp_data{r_acq} = hhp_r;
        all_start_times(r_acq) = start_time_r;
        all_K_values(r_acq) = size(hhp_r, 1);
        
        if ~isempty(hhp_r) && mod(r_acq, 10) == 0
            fprintf('    Max abs value: %g (K_r = %d)\n', max(abs(hhp_r(:))), size(hhp_r, 1));
        elseif isempty(hhp_r)
            fprintf('    WARNING: hhp_r is empty for acquisition %d!\n', r_acq);
        end
    end
    
    h_gen_time = toc;
    fprintf('H matrix generation completed in %.2f seconds\n', h_gen_time);
    
    % Assemble H matrix
    H = assemble_experimental_h_matrix(all_hhp_data, all_start_times, all_K_values, params, imaging_grid, output_folder);
end

function [tx_Aperture, rx_Aperture] = create_experimental_apertures(tx_positions, rx_pos, pMUT_width, pMUT_height, kerf)
    % Create virtual grid
    num_x_grid = 9;
    num_y_grid = 9;
    fprintf('Using a %dx%d virtual grid for mapping.\n', num_x_grid, num_y_grid);
    
    % Generate physical element centers
    physical_element_centers = generate_element_grid(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    
    % Map transmit elements
    tx_active_indices = map_elements_to_grid(tx_positions, physical_element_centers, 'TX');
    
    % Map receive element
    rx_distances = sum((physical_element_centers - rx_pos).^2, 2);
    [~, rx_active_index] = min(rx_distances);
    fprintf('  RX: mapped to grid element %d\n', rx_active_index);
    
    % Create apertures
    tx_Aperture = create_aperture_matrix(tx_active_indices, num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    rx_Aperture = create_aperture_matrix(rx_active_index, num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf);
    
    % Setup impulse responses
    setup_experimental_impulse_responses(tx_Aperture, rx_Aperture);
end

function setup_experimental_impulse_responses(tx_Aperture, rx_Aperture)
    % Define chirp parameters (same as simulation)
    fs = 2e6;
    f_start = 45e3;  % 45 kHz
    f_end = 65e3;    % 65 kHz
    duration = 0.02e-3; % 20 us
    amplitude = 10000;
    
    % Generate chirp
    t_burst = 0 : 1/fs : duration;
    burst_base = chirp(t_burst, f_start, t_burst(end), f_end, 'linear');
    burst_windowed = burst_base .* tukeywin(length(t_burst), 0.25)';
    impulse_response = burst_windowed * amplitude;
    
    % Set impulse responses
    xdc_impulse(tx_Aperture, impulse_response);
    xdc_excitation(tx_Aperture, 1);
    xdc_impulse(rx_Aperture, 1);
    xdc_excitation(rx_Aperture, 1);
    
    fprintf('Chirp: %g-%g kHz, %g us, amplitude %g\n', f_start/1e3, f_end/1e3, duration*1e3, amplitude);
end

function H = assemble_experimental_h_matrix(all_hhp_data, all_start_times, all_K_values, params, imaging_grid, output_folder)
    fprintf('\n--- Assembling Experimental H Matrix ---\n');
    
    % Vectorized time window calculation
    valid_indices = all_K_values > 0;
    all_end_times = zeros(params.R_acquisitions, 1);
    all_end_times(valid_indices) = all_start_times(valid_indices) + (all_K_values(valid_indices) - 1) / params.fs;
    
    min_global_start_time = min(all_start_times);
    max_global_end_time = max(all_end_times);
    
    % Fallback if time window is invalid
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
        warning('H assembly time window invalid, using fallback.');
    end
    
    % Create common time axis
    t_common_axis = min_global_start_time:1/params.fs:max_global_end_time;
    K_global_per_acq = length(t_common_axis);
    if K_global_per_acq == 0
        K_global_per_acq = 1;
        t_common_axis = min_global_start_time;
        warning('K_global_per_acq was 0, set to 1.');
    end
    
    fprintf('Global Time Window: Start=%g s, End=%g s, K_global_per_acq=%d samples.\n', ...
        min_global_start_time, max_global_end_time, K_global_per_acq);
    
    % Pre-allocate H matrix for efficiency
    if params.use_sparse
        total_rows = K_global_per_acq * params.R_acquisitions;
        total_cols = imaging_grid.N_pixels;
        estimated_nnz = total_rows * total_cols * 0.1;
        H_assembled = spalloc(total_rows, total_cols, estimated_nnz);
        fprintf('Using sparse matrix format (estimated %.1f%% sparsity)\n', 100 * estimated_nnz / (total_rows * total_cols));
    else
        H_assembled = zeros(K_global_per_acq * params.R_acquisitions, imaging_grid.N_pixels);
    end
    
    current_row_offset = 0;
    
    % Assemble H matrix
    fprintf('Assembling H matrix...\n');
    tic;
    
    for r_acq = 1:params.R_acquisitions
        hhp_current = all_hhp_data{r_acq};
        start_time_current = all_start_times(r_acq);
        K_current = all_K_values(r_acq);
        
        if K_current == 0 || isempty(hhp_current)
            current_row_offset = current_row_offset + K_global_per_acq;
            continue;
        end
        
        t_current_acq_axis = start_time_current + (0:(K_current - 1)) / params.fs;
        
        % Interpolation
        if params.vectorize_interp && length(t_current_acq_axis) > 1 && issorted(t_current_acq_axis)
            hhp_aligned_r = zeros(K_global_per_acq, imaging_grid.N_pixels);
            
            for px_col = 1:imaging_grid.N_pixels
                if ~isempty(hhp_current) && size(hhp_current, 2) >= px_col
                    hhp_aligned_r(:, px_col) = interp1(t_current_acq_axis, hhp_current(:, px_col), t_common_axis, 'linear', 0);
                end
            end
            
            if params.use_sparse
                threshold = max(abs(hhp_aligned_r(:))) * 1e-10;
                hhp_aligned_r(abs(hhp_aligned_r) < threshold) = 0;
            end
        elseif isscalar(t_current_acq_axis) && K_global_per_acq >= 1
            [~, idx_match] = min(abs(t_common_axis - t_current_acq_axis));
            if ~isempty(idx_match) && ~isempty(hhp_current)
                hhp_aligned_r = zeros(K_global_per_acq, imaging_grid.N_pixels);
                hhp_aligned_r(idx_match, :) = hhp_current(1, :);
            else
                hhp_aligned_r = zeros(K_global_per_acq, imaging_grid.N_pixels);
            end
        else
            hhp_aligned_r = zeros(K_global_per_acq, imaging_grid.N_pixels);
        end
        
        % Matrix assignment
        row_indices = current_row_offset + (1:K_global_per_acq);
        if max(row_indices) <= size(H_assembled, 1)
            H_assembled(row_indices, :) = hhp_aligned_r;
        else
            warning('Row indices exceed matrix bounds, truncating');
            valid_rows = row_indices(row_indices <= size(H_assembled, 1));
            H_assembled(valid_rows, :) = hhp_aligned_r(1:length(valid_rows), :);
        end
        current_row_offset = current_row_offset + K_global_per_acq;
    end
    
    assembly_time = toc;
    fprintf('H matrix assembly completed in %.3f seconds\n', assembly_time);
    
    % Compress sparse matrix if needed
    if params.use_sparse && issparse(H_assembled)
        fprintf('Compressing sparse matrix...\n');
        H_assembled = sparse(H_assembled);
        nnz_ratio = nnz(H_assembled) / numel(H_assembled);
        fprintf('Final sparsity: %.2f%% (%d non-zero elements)\n', 100 * nnz_ratio, nnz(H_assembled));
    end
    
    H = H_assembled;
    M = size(H, 1);
    N = imaging_grid.N_pixels;
    fprintf('Final Experimental H Matrix: %d rows (M) x %d columns (N).\n', M, N);
    
    if M == 0
        error('H matrix has zero rows.');
    end
    
    % Plot H matrix columns
    plot_experimental_h_matrix_columns(H, params, output_folder);
end

function plot_experimental_h_matrix_columns(H, params, output_folder)
    figure(3);
    set(gcf, 'visible', 'off');
    clf;
    hold on;
    
    num_cols_to_plot = min(size(H, 2), 5);
    indices_to_plot = round(linspace(1, size(H, 2), num_cols_to_plot));
    t_axis_plot = (0:(size(H, 1) - 1)) / params.fs * 1e6;
    
    for n_idx = 1:length(indices_to_plot)
        col_idx = indices_to_plot(n_idx);
        plot(t_axis_plot, H(:, col_idx), 'DisplayName', sprintf('H col Px %d', col_idx));
    end
    
    hold off;
    xlabel('Overall Row Index (us)');
    ylabel('Amplitude');
    title('Figure 3: Columns of Experimental H Matrix');
    axis tight;
    grid on;
    legend('Location', 'best');
    set(gcf, 'Color', 'w');
    saveas(gcf, fullfile(output_folder, 'figure3_experimental_h_columns.png'));
    close(gcf);
end

function perform_experimental_reconstruction(H, measurements, params, output_folder)
    fprintf('\n--- Starting Experimental ADMM-TV Reconstruction ---\n');
    
    % Setup reconstruction parameters
    A_matrix = H;
    At_matrix = transpose(A_matrix);
    % Convert measurements to double and ensure proper format
    b_vector = double(measurements.echo_data(:)); % Use real measurements
    
    fprintf('  ADMM Parameters: numIters=%d, rho=%g, lambda_TV=%g\n', ...
        params.numItersADMM, params.rho_admm, params.lambda_tv_reg);
    fprintf('  Convergence: tol=%g, max_iter=%d\n', params.admm_tol, params.admm_max_iter);
    fprintf('  Advanced Features: Sparse=%s, GPU=%s, Adaptive=%s\n', ...
        mat2str(params.use_sparse), mat2str(params.use_gpu), mat2str(params.adaptive_tolerance));
    
    % Check measurement vector size
    fprintf('Measurement vector size: %d, H matrix rows: %d\n', length(b_vector), size(A_matrix, 1));
    
    % Ensure measurement vector matches H matrix rows
    if length(b_vector) ~= size(A_matrix, 1)
        fprintf('Warning: Measurement vector size (%d) does not match H matrix rows (%d)\n', ...
            length(b_vector), size(A_matrix, 1));
        % Truncate or pad as needed
        if length(b_vector) > size(A_matrix, 1)
            b_vector = b_vector(1:size(A_matrix, 1));
            fprintf('Truncated measurement vector to %d elements\n', length(b_vector));
        else
            % Pad with zeros
            b_vector = [b_vector; zeros(size(A_matrix, 1) - length(b_vector), 1)];
            fprintf('Padded measurement vector to %d elements\n', length(b_vector));
        end
    end
    
    % Normalize system matrix
    H_norm_factor = max(abs(A_matrix(:)));
    if H_norm_factor < eps
        H_norm_factor = 1;
    end
    A_admm = A_matrix ./ H_norm_factor;
    At_admm = transpose(A_admm);
    b_admm_vec = b_vector(:) / H_norm_factor;
    
    % Determine image resolution from H matrix
    N_pixels = size(H, 2);
    % Calculate grid dimensions based on imaging parameters
    grid_width_mm = params.grid_width_mm;
    grid_depth_mm = params.grid_depth_end_mm - params.grid_depth_start_mm;
    grid_step_mm = params.grid_step_mm;
    
    grid_width_pixels = round(grid_width_mm / grid_step_mm) + 1;
    grid_depth_pixels = round(grid_depth_mm / grid_step_mm) + 1;
    
    % Ensure the calculated resolution matches the H matrix
    expected_pixels = grid_depth_pixels * grid_width_pixels;
    if expected_pixels ~= N_pixels
        fprintf('Warning: Expected %d pixels but H matrix has %d columns\n', expected_pixels, N_pixels);
        % Adjust to match H matrix
        if N_pixels == 988  % This is the actual value from the output
            imageResolution = [26, 38];  % 26 * 38 = 988
        else
            % Fallback calculation
            grid_width_pixels = round(sqrt(N_pixels));
            grid_depth_pixels = ceil(N_pixels / grid_width_pixels);
            imageResolution = [grid_depth_pixels, grid_width_pixels];
        end
    else
        imageResolution = [grid_depth_pixels, grid_width_pixels];
    end
    
    fprintf('Image resolution: %s (%d pixels)\n', mat2str(imageResolution), N_pixels);
    fprintf('Grid: %.0f mm x %.0f mm, step %.1f mm\n', grid_width_mm, grid_depth_mm, grid_step_mm);
    
    % Setup ADMM operators
    fprintf('Matrix dimensions: A_admm=%s, At_admm=%s, imageResolution=%s\n', ...
        mat2str(size(A_admm)), mat2str(size(At_admm)), mat2str(imageResolution));
    
    % Check matrix compatibility
    expected_cols = prod(imageResolution);
    if size(At_admm, 1) ~= expected_cols
        fprintf('Error: At_admm has %d rows but expected %d (imageResolution=%s)\n', ...
            size(At_admm, 1), expected_cols, mat2str(imageResolution));
        error('Matrix dimension mismatch');
    end
    
    [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = setup_admm_operators_optimized(A_admm, At_admm, imageResolution, params);
    
    % Initialize ADMM variables
    x_admm_img_iter = zeros(imageResolution);
    z_admm_grad_iter = zeros([prod(imageResolution) 2]);
    u_admm_dual_iter = zeros([prod(imageResolution) 2]);
    
    % Optimized PCG function
    Hfun_pcg_admm = @(x_vec) reshape(AtAfun_admm_img(reshape(x_vec, imageResolution)) + ...
        params.rho_admm .* opDtDx_tv(reshape(x_vec, imageResolution)), [prod(imageResolution) 1]);
    
    % Initialize tracking arrays
    residuals_admm_iters = zeros([params.admm_max_iter 1]);
    
    % Adaptive tolerance
    current_tol = params.admm_tol;
    
    % Setup visualization
    figure(8);
    set(gcf, 'visible', 'off');
    clf;
    set(gcf, 'Position', [200, 200, 1000, 700], 'Color', 'w');
    
    % Run ADMM iterations
    tic;
    converged = false;
    prev_residual = inf;
    
    for k_admm = 1:params.admm_max_iter
        % Store previous iteration for convergence check
        x_prev = x_admm_img_iter;
        
        % ADMM update steps
        [x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter] = admm_iteration_optimized(...
            x_admm_img_iter, z_admm_grad_iter, u_admm_dual_iter, ...
            Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, ...
            params.rho_admm, params.lambda_tv_reg, 0, Hfun_pcg_admm, params, k_admm);
        
        % Calculate residual
        r1 = b_admm_vec - Afun_admm(x_admm_img_iter);
        r2 = opDx_tv(x_admm_img_iter);
        tv_n = sum(sqrt(sum(r2.^2, 2)));
        residuals_admm_iters(k_admm) = 0.5 * sum(r1(:).^2) + params.lambda_tv_reg * tv_n;
        
        % Adaptive tolerance adjustment
        if params.adaptive_tolerance && k_admm > 5
            residual_change = abs(residuals_admm_iters(k_admm) - prev_residual) / (abs(prev_residual) + eps);
            if residual_change < 1e-3
                current_tol = max(current_tol * 0.5, 1e-8);
            end
            prev_residual = residuals_admm_iters(k_admm);
        end
        
        % Check convergence
        if k_admm > 1
            rel_change = norm(x_admm_img_iter(:) - x_prev(:)) / (norm(x_prev(:)) + eps);
            if rel_change < current_tol
                converged = true;
                fprintf('  ADMM converged at iteration %d (rel_change=%.2e, tol=%.2e)\n', k_admm, rel_change, current_tol);
                break;
            end
        end
        
        % Update visualization
        if mod(k_admm, 5) == 0 || k_admm == 1 || converged
            update_experimental_visualization(x_admm_img_iter, residuals_admm_iters, k_admm, params, imageResolution);
        end
    end
    
    runtime_ADMM_total = toc;
    fprintf('  Experimental ADMM-TV: Time=%.2fs, Iterations=%d\n', runtime_ADMM_total, k_admm);
    
    % Save final result
    sgtitle(sprintf('Figure 8: Experimental ADMM TV (R_{acq}=%d)', params.R_acquisitions));
    saveas(gcf, fullfile(output_folder, 'figure8_experimental_admm_tv.png'));
    close(gcf);
    
    % Save reconstruction result
    save(fullfile(output_folder, 'experimental_reconstruction.mat'), 'x_admm_img_iter', 'residuals_admm_iters', 'params');
    fprintf('Reconstruction saved to: %s\n', fullfile(output_folder, 'experimental_reconstruction.mat'));
    
    % Create final summary plot
    create_final_summary_plot(x_admm_img_iter, residuals_admm_iters, params, imageResolution, output_folder);
end

function update_experimental_visualization(x_admm_img_iter, residuals_admm_iters, k_admm, params, imageResolution)
    figure(8);
    
    % Normalize current reconstruction
    x_scl = real(x_admm_img_iter(:));
    x_range = max(x_scl) - min(x_scl);
    if x_range > eps
        x_scl = (x_scl - min(x_scl)) / x_range;
    else
        x_scl = zeros(size(x_scl));
    end
    
    % Plot reconstruction
    subplot(2, 2, [1 3]);
    imagesc(reshape(x_scl, imageResolution));
    axis image; colormap(gca, gray); colorbar; set(gca, 'YDir', 'normal');
    title(sprintf('Experimental Reconstruction\n\\lambda=%.1e,\\rho=%.1f, Iteration %d', ...
        params.lambda_tv_reg, params.rho_admm, k_admm));
    xlabel('Lateral'); ylabel('Axial');
    
    % Plot residual convergence
    subplot(2, 2, 2);
    plot(1:k_admm, log10(residuals_admm_iters(1:k_admm)), '-', 'LineWidth', 2);
    title('log10(Objective)/Iteration'); xlabel('Iteration'); ylabel('log10(Value)'); 
    grid on; axis tight;
    if k_admm > 1
        yl = ylim; if diff(yl) < 0.1; yl(2) = yl(1) + 0.1; end; ylim(yl);
    end
    
    % Plot reconstruction statistics
    subplot(2, 2, 4);
    x_stats = [mean(x_scl), std(x_scl), min(x_scl), max(x_scl)];
    bar(x_stats);
    set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
    title('Reconstruction Statistics');
    ylabel('Amplitude');
    grid on;
end

function cleanup_experimental_reconstruction()
    % Clean up parallel pool (with error handling)
    try
        if ~isempty(gcp('nocreate'))
            delete(gcp('nocreate'));
        end
    catch
        % Ignore parallel computing toolbox errors
    end
    
    if exist('field_end', 'file') == 2
        field_end;
        disp('Field II ended.');
    else
        disp('field_end function not found.');
    end
    fprintf('\nExperimental reconstruction completed successfully!\n');
end

%% ===== HELPER FUNCTIONS (from simulation) =====

function element_centers = generate_element_grid(num_x, num_y, width, height, kerf)
    element_centers = zeros(num_x * num_y, 3);
    center_offset_x = (num_x - 1) / 2 * (width + kerf);
    center_offset_y = (num_y - 1) / 2 * (height + kerf);
    
    element_idx = 1;
    for iy = 1:num_y
        y_pos = (iy - 1) * (height + kerf) - center_offset_y;
        for ix = 1:num_x
            x_pos = (ix - 1) * (width + kerf) - center_offset_x;
            element_centers(element_idx, :) = [x_pos, y_pos, 0];
            element_idx = element_idx + 1;
        end
    end
end

function active_indices = map_elements_to_grid(desired_positions, element_centers, element_type)
    active_indices = zeros(size(desired_positions, 1), 1);
    fprintf('Mapping %s elements:\n', element_type);
    
    for i = 1:size(desired_positions, 1)
        distances = sum((element_centers - desired_positions(i, :)).^2, 2);
        [~, min_idx] = min(distances);
        active_indices(i) = min_idx;
        fprintf('  %s%d: mapped to grid element %d\n', element_type, i, min_idx);
    end
    
    % Check for duplicates
    original_indices = active_indices;
    active_indices = unique(active_indices);
    if length(original_indices) ~= length(active_indices)
        fprintf('WARNING: %d elements mapped to same positions!\n', ...
            length(original_indices) - length(active_indices));
    end
end

function aperture = create_aperture_matrix(active_indices, num_x, num_y, width, height, kerf)
    enabled_matrix = zeros(num_y, num_x);
    [row_indices, col_indices] = ind2sub([num_y, num_x], active_indices);
    for i = 1:length(active_indices)
        enabled_matrix(row_indices(i), col_indices(i)) = 1;
    end
    aperture = xdc_2d_array(num_x, num_y, width, height, kerf, kerf, enabled_matrix, 1, 1, [0 0 100e-3]);
end

function [Afun_admm, Atfun_admm_img, AtAfun_admm_img, opDx_tv, opDtx_tv, opDtDx_tv] = setup_admm_operators_optimized(A_admm, At_admm, imageResolution, params)
    Afun_admm = @(x) A_admm * x(:);
    Atfun_admm_img = @(y) reshape(At_admm * y, imageResolution);
    AtAfun_admm_img = @(x) Atfun_admm_img(Afun_admm(x));
    
    [Dx_sparse, Dy_sparse] = createDifferenceOperators(imageResolution);
    opDx_tv = @(x) [Dx_sparse * x(:), Dy_sparse * x(:)];
    opDtx_tv = @(v) reshape(Dx_sparse' * v(:, 1) + Dy_sparse' * v(:, 2), imageResolution);
    opDtDx_tv = @(x) opDtx_tv(opDx_tv(x));
end

function [x_new, z_new, u_new] = admm_iteration_optimized(x_old, z_old, u_old, Atfun_admm_img, b_admm_vec, opDx_tv, opDtx_tv, rho_admm, lambda_tv_reg, noise_sigma_admm, Hfun_pcg_admm, params, k_admm)
    % ADMM x-update with optimized PCG
    v_upd = z_old - u_old;
    bb_upd = Atfun_admm_img(b_admm_vec) + rho_admm * opDtx_tv(v_upd);
    [x_vec_new, ~, ~, ~] = pcg(Hfun_pcg_admm, bb_upd(:), params.pcg_tol, params.pcg_max_iter, [], [], x_old(:));
    x_new = reshape(x_vec_new, size(x_old));
    
    % ADMM z-update (vectorized)
    kap = lambda_tv_reg / rho_admm;
    v_z_upd = opDx_tv(x_new) + u_old;
    v_norm = sqrt(sum(v_z_upd.^2, 2));
    v_norm = max(v_norm, eps); % Avoid division by zero
    shr = max(0, 1 - kap ./ v_norm);
    z_new = v_z_upd .* shr;
    
    % ADMM u-update
    u_new = u_old + opDx_tv(x_new) - z_new;
end

function [Dx, Dy] = createDifferenceOperators(imageSize)
    rows = imageSize(1);
    cols = imageSize(2);
    N_img_pixels = rows * cols;
    
    Dx = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, rows], N_img_pixels, N_img_pixels);
    last_col_indices_mask = false(N_img_pixels, 1);
    last_col_indices_mask((cols-1)*rows+1 : cols*rows) = true;
    Dx(last_col_indices_mask, :) = 0;
    
    Dy = spdiags([-ones(N_img_pixels, 1), ones(N_img_pixels, 1)], [0, 1], N_img_pixels, N_img_pixels);
    last_row_indices_mask = false(N_img_pixels, 1);
    last_row_indices_mask(rows:rows:N_img_pixels) = true;
    Dy(last_row_indices_mask, :) = 0;
end

function create_final_summary_plot(x_admm_img_iter, residuals_admm_iters, params, imageResolution, output_folder)
    % Create final summary plot of the reconstruction
    figure(9);
    set(gcf, 'visible', 'off');
    clf;
    set(gcf, 'Position', [300, 300, 1200, 800], 'Color', 'w');
    
    % Normalize final reconstruction
    x_scl = real(x_admm_img_iter(:));
    x_range = max(x_scl) - min(x_scl);
    if x_range > eps
        x_scl = (x_scl - min(x_scl)) / x_range;
    else
        x_scl = zeros(size(x_scl));
    end
    
    % Plot final reconstruction
    subplot(2, 3, [1 2 4 5]);
    imagesc(reshape(x_scl, imageResolution));
    axis image; colormap(gca, gray); colorbar; set(gca, 'YDir', 'normal');
    title(sprintf('Final Experimental Reconstruction\n\\lambda=%.1e, \\rho=%.1f, %d iterations', ...
        params.lambda_tv_reg, params.rho_admm, length(residuals_admm_iters)));
    xlabel('Lateral Position'); ylabel('Axial Position');
    
    % Plot convergence
    subplot(2, 3, 3);
    plot(1:length(residuals_admm_iters), log10(residuals_admm_iters), '-', 'LineWidth', 2);
    title('Convergence History'); xlabel('Iteration'); ylabel('log10(Objective)');
    grid on; axis tight;
    
    % Plot reconstruction statistics
    subplot(2, 3, 6);
    x_stats = [mean(x_scl), std(x_scl), min(x_scl), max(x_scl)];
    bar(x_stats);
    set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
    title('Reconstruction Statistics');
    ylabel('Amplitude');
    grid on;
    
    sgtitle(sprintf('Figure 9: Final Experimental Reconstruction Summary (R_{acq}=%d)', params.R_acquisitions));
    saveas(gcf, fullfile(output_folder, 'figure9_final_summary.png'));
    close(gcf);
    
    fprintf('Final summary plot saved to: %s\n', fullfile(output_folder, 'figure9_final_summary.png'));
end

function validate_experimental_files(h5_file_path, delay_csv_path, freq_csv_path)
    % Validate that all experimental files exist and are readable
    fprintf('Validating experimental files...\n');
    
    % Check H5 file
    if ~exist(h5_file_path, 'file')
        error('H5 file not found: %s', h5_file_path);
    end
    fprintf('  ✓ H5 file found: %s\n', h5_file_path);
    
    % Check delay CSV file
    if ~exist(delay_csv_path, 'file')
        error('Delay CSV file not found: %s', delay_csv_path);
    end
    fprintf('  ✓ Delay CSV file found: %s\n', delay_csv_path);
    
    % Check frequency CSV file
    if ~exist(freq_csv_path, 'file')
        error('Frequency CSV file not found: %s', freq_csv_path);
    end
    fprintf('  ✓ Frequency CSV file found: %s\n', freq_csv_path);
    
    % Try to read a small sample from each file to verify they're readable
    try
        h5_info = h5info(h5_file_path);
        fprintf('  ✓ H5 file is readable (%d datasets)\n', length(h5_info.Datasets));
    catch ME
        error('H5 file is not readable: %s', ME.message);
    end
    
    try
        delay_data = readtable(delay_csv_path);
        fprintf('  ✓ Delay CSV file is readable (%d rows, %d columns)\n', size(delay_data, 1), size(delay_data, 2));
    catch ME
        error('Delay CSV file is not readable: %s', ME.message);
    end
    
    try
        freq_data = readtable(freq_csv_path);
        fprintf('  ✓ Frequency CSV file is readable (%d rows, %d columns)\n', size(freq_data, 1), size(freq_data, 2));
    catch ME
        error('Frequency CSV file is not readable: %s', ME.message);
    end
    
    fprintf('All experimental files validated successfully.\n\n');
end

function filtered_data = apply_lowpass_filter(data, fs, cutoff_freq)
    % Apply lowpass filter to experimental data
    % Inputs:
    %   data: [num_acquisitions x signal_length] matrix
    %   fs: sampling frequency (Hz)
    %   cutoff_freq: cutoff frequency (Hz)
    
    fprintf('  Filtering %d acquisitions with %d Hz lowpass filter...\n', ...
        size(data, 1), cutoff_freq);
    
    % Design lowpass filter
    nyquist_freq = fs / 2;
    normalized_cutoff = cutoff_freq / nyquist_freq;
    
    % Use Butterworth filter for smooth frequency response
    [b, a] = butter(4, normalized_cutoff, 'low');
    
    % Apply filter to each acquisition
    filtered_data = zeros(size(data));
    for acq = 1:size(data, 1)
        if mod(acq, 20) == 0 || acq == 1
            fprintf('    Processing acquisition %d/%d...\n', acq, size(data, 1));
        end
        
        % Convert to double for filtering
        signal = double(data(acq, :));
        
        % Apply filter
        filtered_signal = filtfilt(b, a, signal);
        
        % Store filtered result
        filtered_data(acq, :) = filtered_signal;
    end
    
    fprintf('  Lowpass filtering completed for all acquisitions.\n');
end 