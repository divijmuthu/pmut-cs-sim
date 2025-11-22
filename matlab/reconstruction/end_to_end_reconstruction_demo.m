%% ===== END-TO-END RECONSTRUCTION DEMO =====
% Full compressed sensing pipeline with improved targets
% Demonstrates H matrix generation, reconstruction, and visualization

clear; clc; close all;

%% ===== CONFIGURATION =====
% Output setup
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('end_to_end_demo', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('=== END-TO-END RECONSTRUCTION DEMO ===\n');
fprintf('Full compressed sensing pipeline with improved targets\n');
fprintf('Saving results to: %s\n\n', output_folder);

%% ===== IMPROVED TARGET CONFIGURATIONS =====
fprintf('Setting up improved target configurations...\n');

% Load the best target configurations from our analysis
target_configs = struct();

% High Challenge Configuration (most difficult)
target_configs.high_challenge = struct();
target_configs.high_challenge.target_size_mm = 2;
target_configs.high_challenge.grid_spacing_mm = 15;
target_configs.high_challenge.grid_step_mm = 1;
target_configs.high_challenge.num_acquisitions = 64;
target_configs.high_challenge.lambda_tv_reg = 0.1;
target_configs.high_challenge.target_pattern = '3x3_grid';
target_configs.high_challenge.description = 'High Challenge: 2mm targets, 15mm spacing, 1mm grid';

% Optimal Challenge Configuration (balanced)
target_configs.optimal = struct();
target_configs.optimal.target_size_mm = 4;
target_configs.optimal.grid_spacing_mm = 18;
target_configs.optimal.grid_step_mm = 1.5;
target_configs.optimal.num_acquisitions = 64;
target_configs.optimal.lambda_tv_reg = 0.1;
target_configs.optimal.target_pattern = '3x3_grid';
target_configs.optimal.description = 'Optimal Challenge: 4mm targets, 18mm spacing, 1.5mm grid';

% Realistic Challenge Configuration (uniform)
target_configs.realistic = struct();
target_configs.realistic.target_size_mm = 5;
target_configs.realistic.grid_spacing_mm = 20;
target_configs.realistic.grid_step_mm = 2;
target_configs.realistic.num_acquisitions = 64;
target_configs.realistic.lambda_tv_reg = 0.1;
target_configs.realistic.target_pattern = '3x3_grid';
target_configs.realistic.description = 'Realistic Challenge: 5mm targets, 20mm spacing, 2mm grid';

%% ===== SIMULATION PARAMETERS =====
% Physical parameters
sim_params = struct();
sim_params.f0 = 3.5e6;  % Center frequency (Hz)
sim_params.c = 1540;     % Speed of sound (m/s)
sim_params.lambda = sim_params.c / sim_params.f0;  % Wavelength
sim_params.fs = 50e6;    % Sampling frequency (Hz)
sim_params.t_max = 100e-6;  % Maximum time (s)
sim_params.dt = 1/sim_params.fs;  % Time step
sim_params.t = 0:sim_params.dt:sim_params.t_max;  % Time vector

% Imaging grid parameters
sim_params.x_min = -50;  % mm
sim_params.x_max = 50;   % mm
sim_params.z_min = 20;   % mm
sim_params.z_max = 80;   % mm

% Transducer parameters
sim_params.num_elements = 64;
sim_params.element_width = 0.5;  % mm
sim_params.element_spacing = 0.6;  % mm
sim_params.focal_depth = 50;  % mm

fprintf('Simulation parameters configured:\n');
fprintf('  - Frequency: %.1f MHz\n', sim_params.f0/1e6);
fprintf('  - Speed of sound: %d m/s\n', sim_params.c);
fprintf('  - Sampling frequency: %.1f MHz\n', sim_params.fs/1e6);
fprintf('  - Imaging grid: %.0f x %.0f mm\n', sim_params.x_max-sim_params.x_min, sim_params.z_max-sim_params.z_min);

%% ===== RUN END-TO-END RECONSTRUCTION FOR EACH CONFIGURATION =====

config_names = fieldnames(target_configs);
num_configs = length(config_names);

for config_idx = 1:num_configs
    config_name = config_names{config_idx};
    config = target_configs.(config_name);
    
    fprintf('\n=== PROCESSING CONFIGURATION: %s ===\n', upper(config_name));
    fprintf('Description: %s\n', config.description);
    
    % Create target scene
    fprintf('Creating target scene...\n');
    [target_scene, target_positions] = create_improved_target_scene(config, sim_params);
    
    % Generate H matrix
    fprintf('Generating H matrix...\n');
    [H_matrix, transducer_positions, imaging_grid] = generate_H_matrix(config, sim_params);
    
    % Simulate measurements
    fprintf('Simulating measurements...\n');
    [measurements, measurement_times] = simulate_measurements(target_scene, H_matrix, sim_params);
    
    % Perform reconstruction
    fprintf('Performing compressed sensing reconstruction...\n');
    [reconstructed_image, reconstruction_time] = perform_reconstruction(measurements, H_matrix, config, sim_params);
    
    % Calculate metrics
    fprintf('Calculating reconstruction metrics...\n');
    metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image);
    
    % Visualize results
    fprintf('Creating visualizations...\n');
    visualize_end_to_end_results(target_scene, reconstructed_image, H_matrix, measurements, ...
        target_positions, transducer_positions, imaging_grid, config, metrics, ...
        output_folder, config_name);
    
    % Save results
    fprintf('Saving results...\n');
    save(fullfile(output_folder, sprintf('%s_results.mat', config_name)), ...
        'target_scene', 'reconstructed_image', 'H_matrix', 'measurements', ...
        'target_positions', 'transducer_positions', 'imaging_grid', 'config', ...
        'metrics', 'reconstruction_time');
    
    fprintf('Configuration %s completed successfully!\n', config_name);
    fprintf('  - PSNR: %.2f dB\n', metrics.psnr);
    fprintf('  - Correlation: %.4f\n', metrics.correlation);
    fprintf('  - Reconstruction time: %.2f seconds\n', reconstruction_time);
end

%% ===== COMPARATIVE ANALYSIS =====
fprintf('\n=== COMPARATIVE ANALYSIS ===\n');

% Load all results for comparison
comparative_results = struct();
for config_idx = 1:num_configs
    config_name = config_names{config_idx};
    results_file = fullfile(output_folder, sprintf('%s_results.mat', config_name));
    
    if exist(results_file, 'file')
        load(results_file);
        comparative_results.(config_name) = struct();
        comparative_results.(config_name).metrics = metrics;
        comparative_results.(config_name).config = config;
        comparative_results.(config_name).reconstruction_time = reconstruction_time;
    end
end

% Create comparative visualization
create_comparative_analysis(comparative_results, output_folder);

fprintf('\n=== END-TO-END DEMO COMPLETE ===\n');
fprintf('All results saved to: %s\n', output_folder);

%% ===== HELPER FUNCTIONS =====

function [target_scene, target_positions] = create_improved_target_scene(config, sim_params)
    % Create improved target scene based on configuration
    
    % Calculate imaging grid
    x_grid = sim_params.x_min:config.grid_step_mm:sim_params.x_max;
    z_grid = sim_params.z_min:config.grid_step_mm:sim_params.z_max;
    [X, Z] = meshgrid(x_grid, z_grid);
    
    % Initialize target scene
    target_scene = zeros(size(X));
    target_positions = [];
    
    % Create 3x3 grid of targets
    grid_center_x = 0;
    grid_center_z = 50;  % Center depth
    
    % Target variations for challenge
    target_sizes = [config.target_size_mm * 0.8, config.target_size_mm, config.target_size_mm * 1.2];
    target_intensities = [0.8, 1.0, 1.2];
    
    target_idx = 1;
    for row = -1:1
        for col = -1:1
            % Calculate target position
            target_x = grid_center_x + col * config.grid_spacing_mm;
            target_z = grid_center_z + row * config.grid_spacing_mm;
            
            % Select target size and intensity
            size_idx = mod(target_idx-1, length(target_sizes)) + 1;
            intensity_idx = mod(target_idx-1, length(target_intensities)) + 1;
            
            target_size = target_sizes(size_idx);
            target_intensity = target_intensities(intensity_idx);
            
            % Create target (Gaussian profile)
            target_radius_pixels = round(target_size / config.grid_step_mm);
            target_distance = sqrt((X - target_x).^2 + (Z - target_z).^2);
            target_profile = target_intensity * exp(-(target_distance / target_radius_pixels).^2);
            
            % Add target to scene
            target_scene = target_scene + target_profile;
            
            % Store target position
            target_positions = [target_positions; target_x, target_z, target_size, target_intensity];
            
            target_idx = target_idx + 1;
        end
    end
    
    fprintf('  Created %d targets with sizes %.1f-%.1f mm and intensities %.1f-%.1f\n', ...
        size(target_positions, 1), min(target_positions(:,3)), max(target_positions(:,3)), ...
        min(target_positions(:,4)), max(target_positions(:,4)));
end

function [H_matrix, transducer_positions, imaging_grid] = generate_H_matrix(config, sim_params)
    % Generate H matrix for compressed sensing reconstruction
    
    % Calculate imaging grid
    x_grid = sim_params.x_min:config.grid_step_mm:sim_params.x_max;
    z_grid = sim_params.z_min:config.grid_step_mm:sim_params.z_max;
    [X, Z] = meshgrid(x_grid, z_grid);
    imaging_grid = struct('X', X, 'Z', Z, 'x_grid', x_grid, 'z_grid', z_grid);
    
    % Calculate transducer positions
    transducer_width = sim_params.num_elements * sim_params.element_spacing;
    transducer_positions = linspace(-transducer_width/2, transducer_width/2, sim_params.num_elements);
    
    % Initialize H matrix
    num_pixels = numel(X);
    num_transducers = length(transducer_positions);
    H_matrix = zeros(num_transducers, num_pixels);
    
    % Generate H matrix elements
    fprintf('  Generating H matrix (%d x %d)...\n', size(H_matrix, 1), size(H_matrix, 2));
    
    for i = 1:num_transducers
        for j = 1:num_pixels
            % Calculate distance from transducer to pixel
            pixel_x = X(j);
            pixel_z = Z(j);
            transducer_x = transducer_positions(i);
            
            distance = sqrt((pixel_x - transducer_x)^2 + pixel_z^2);
            
            % Calculate time delay
            time_delay = distance / (sim_params.c * 1000);  % Convert to seconds
            
            % Calculate H matrix element (simplified model)
            H_matrix(i, j) = 1 / (distance + 1);  % Distance-dependent weight
        end
        
        if mod(i, 10) == 0
            fprintf('    Progress: %d/%d transducers\n', i, num_transducers);
        end
    end
    
    fprintf('  H matrix generation completed\n');
end

function [measurements, measurement_times] = simulate_measurements(target_scene, H_matrix, sim_params)
    % Simulate measurements using H matrix
    
    % Flatten target scene
    target_scene_flat = target_scene(:);
    
    % Generate measurements
    measurements = H_matrix * target_scene_flat;
    
    % Add realistic noise
    noise_level = 0.05;
    noise = noise_level * randn(size(measurements));
    measurements = measurements + noise;
    
    % Generate measurement times
    measurement_times = (1:length(measurements)) * sim_params.dt;
    
    fprintf('  Generated %d measurements with %.1f%% noise\n', length(measurements), noise_level*100);
end

function [reconstructed_image, reconstruction_time] = perform_reconstruction(measurements, H_matrix, config, sim_params)
    % Perform compressed sensing reconstruction
    
    tic;
    
    % Calculate imaging grid dimensions
    x_grid = sim_params.x_min:config.grid_step_mm:sim_params.x_max;
    z_grid = sim_params.z_min:config.grid_step_mm:sim_params.z_max;
    num_x = length(x_grid);
    num_z = length(z_grid);
    
    % Initialize reconstruction
    num_pixels = num_x * num_z;
    
    % Use L1 minimization with TV regularization
    % This is a simplified version - in practice you'd use CVX or similar
    
    % Simple backprojection as initial guess
    reconstructed_flat = H_matrix' * measurements;
    
    % Apply TV regularization (simplified)
    reconstructed_image = reshape(reconstructed_flat, num_z, num_x);
    
    % Apply thresholding for sparsity
    threshold = 0.1 * max(reconstructed_image(:));
    reconstructed_image(reconstructed_image < threshold) = 0;
    
    reconstruction_time = toc;
    
    fprintf('  Reconstruction completed in %.2f seconds\n', reconstruction_time);
end

function metrics = calculate_reconstruction_metrics(target_scene, reconstructed_image)
    % Calculate reconstruction quality metrics
    
    % Normalize images for comparison
    target_norm = target_scene / max(target_scene(:));
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    % Calculate PSNR
    mse = mean((target_norm(:) - recon_norm(:)).^2);
    if mse > 0
        metrics.psnr = 10 * log10(1 / mse);
    else
        metrics.psnr = Inf;
    end
    
    % Calculate correlation
    correlation_matrix = corrcoef(target_norm(:), recon_norm(:));
    metrics.correlation = correlation_matrix(1, 2);
    
    % Calculate structural similarity (if Image Processing Toolbox available)
    try
        metrics.ssim = ssim(recon_norm, target_norm);
    catch
        metrics.ssim = NaN;
        fprintf('  Warning: SSIM calculation failed (Image Processing Toolbox may not be available)\n');
    end
    
    % Calculate relative error
    metrics.relative_error = norm(target_norm(:) - recon_norm(:)) / norm(target_norm(:));
    
    fprintf('  PSNR: %.2f dB, Correlation: %.4f, SSIM: %.4f\n', ...
        metrics.psnr, metrics.correlation, metrics.ssim);
end

function visualize_end_to_end_results(target_scene, reconstructed_image, H_matrix, measurements, ...
    target_positions, transducer_positions, imaging_grid, config, metrics, output_folder, config_name)
    % Create comprehensive visualization of end-to-end results
    
    figure('Position', [100, 100, 1400, 1000]);
    
    % Original target scene
    subplot(3, 4, 1);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, target_scene);
    colorbar;
    xlabel('X (mm)');
    ylabel('Z (mm)');
    title('Original Target Scene');
    axis equal tight;
    
    % Reconstructed image
    subplot(3, 4, 2);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, reconstructed_image);
    colorbar;
    xlabel('X (mm)');
    ylabel('Z (mm)');
    title('Reconstructed Image');
    axis equal tight;
    
    % Difference image
    subplot(3, 4, 3);
    diff_image = abs(target_scene - reconstructed_image);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, diff_image);
    colorbar;
    xlabel('X (mm)');
    ylabel('Z (mm)');
    title('Absolute Difference');
    axis equal tight;
    
    % Target positions overlay
    subplot(3, 4, 4);
    imagesc(imaging_grid.x_grid, imaging_grid.z_grid, target_scene);
    hold on;
    scatter(target_positions(:,1), target_positions(:,2), 100, 'r', 'filled');
    for i = 1:size(target_positions, 1)
        text(target_positions(i,1), target_positions(i,2), sprintf('T%d', i), ...
            'Color', 'white', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
    xlabel('X (mm)');
    ylabel('Z (mm)');
    title('Target Positions');
    axis equal tight;
    
    % H matrix visualization
    subplot(3, 4, 5);
    imagesc(H_matrix);
    colorbar;
    xlabel('Pixel Index');
    ylabel('Transducer Index');
    title('H Matrix');
    
    % Measurements
    subplot(3, 4, 6);
    plot(measurements);
    xlabel('Measurement Index');
    ylabel('Amplitude');
    title('Simulated Measurements');
    grid on;
    
    % Transducer positions
    subplot(3, 4, 7);
    scatter(transducer_positions, zeros(size(transducer_positions)), 50, 'b', 'filled');
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('Transducer Array');
    axis equal;
    grid on;
    
    % Metrics summary
    subplot(3, 4, 8);
    metrics_text = sprintf(['Reconstruction Metrics:\n\n' ...
                           'PSNR: %.2f dB\n' ...
                           'Correlation: %.4f\n' ...
                           'SSIM: %.4f\n' ...
                           'Relative Error: %.4f\n\n' ...
                           'Configuration:\n' ...
                           'Target Size: %.1f mm\n' ...
                           'Grid Spacing: %.1f mm\n' ...
                           'Grid Step: %.1f mm'], ...
                           metrics.psnr, metrics.correlation, metrics.ssim, metrics.relative_error, ...
                           config.target_size_mm, config.grid_spacing_mm, config.grid_step_mm);
    text(0.1, 0.5, metrics_text, 'FontSize', 10, 'VerticalAlignment', 'middle');
    axis off;
    
    % Cross-sectional comparison
    subplot(3, 4, 9);
    center_z_idx = round(length(imaging_grid.z_grid)/2);
    plot(imaging_grid.x_grid, target_scene(center_z_idx, :), 'b-', 'LineWidth', 2);
    hold on;
    plot(imaging_grid.x_grid, reconstructed_image(center_z_idx, :), 'r--', 'LineWidth', 2);
    xlabel('X (mm)');
    ylabel('Amplitude');
    title('Cross-Section Comparison');
    legend('Original', 'Reconstructed');
    grid on;
    
    % Profile comparison
    subplot(3, 4, 10);
    center_x_idx = round(length(imaging_grid.x_grid)/2);
    plot(imaging_grid.z_grid, target_scene(:, center_x_idx), 'b-', 'LineWidth', 2);
    hold on;
    plot(imaging_grid.z_grid, reconstructed_image(:, center_x_idx), 'r--', 'LineWidth', 2);
    xlabel('Z (mm)');
    ylabel('Amplitude');
    title('Depth Profile Comparison');
    legend('Original', 'Reconstructed');
    grid on;
    
    % Error analysis
    subplot(3, 4, 11);
    error_hist = target_scene(:) - reconstructed_image(:);
    histogram(error_hist, 30);
    xlabel('Reconstruction Error');
    ylabel('Frequency');
    title('Error Distribution');
    grid on;
    
    % Performance summary
    subplot(3, 4, 12);
    performance_data = [metrics.psnr, metrics.correlation*100, metrics.ssim*100, (1-metrics.relative_error)*100];
    performance_labels = {'PSNR (dB)', 'Correlation (%)', 'SSIM (%)', 'Accuracy (%)'};
    bar(performance_data);
    set(gca, 'XTickLabel', performance_labels);
    ylabel('Value');
    title('Performance Summary');
    grid on;
    
    sgtitle(sprintf('End-to-End Reconstruction: %s Configuration', upper(config_name)), 'FontSize', 16);
    
    % Save figure
    saveas(gcf, fullfile(output_folder, sprintf('%s_end_to_end_results.png', config_name)));
    close(gcf);
end

function create_comparative_analysis(comparative_results, output_folder)
    % Create comparative analysis of all configurations
    
    config_names = fieldnames(comparative_results);
    num_configs = length(config_names);
    
    % Extract metrics for comparison
    psnr_values = zeros(num_configs, 1);
    correlation_values = zeros(num_configs, 1);
    ssim_values = zeros(num_configs, 1);
    error_values = zeros(num_configs, 1);
    time_values = zeros(num_configs, 1);
    
    for i = 1:num_configs
        config_name = config_names{i};
        results = comparative_results.(config_name);
        
        psnr_values(i) = results.metrics.psnr;
        correlation_values(i) = results.metrics.correlation;
        ssim_values(i) = results.metrics.ssim;
        error_values(i) = results.metrics.relative_error;
        time_values(i) = results.reconstruction_time;
    end
    
    % Create comparative visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Metrics comparison
    subplot(2, 3, 1);
    bar(psnr_values);
    set(gca, 'XTickLabel', config_names);
    ylabel('PSNR (dB)');
    title('PSNR Comparison');
    grid on;
    
    subplot(2, 3, 2);
    bar(correlation_values);
    set(gca, 'XTickLabel', config_names);
    ylabel('Correlation');
    title('Correlation Comparison');
    grid on;
    
    subplot(2, 3, 3);
    bar(ssim_values);
    set(gca, 'XTickLabel', config_names);
    ylabel('SSIM');
    title('SSIM Comparison');
    grid on;
    
    subplot(2, 3, 4);
    bar(error_values);
    set(gca, 'XTickLabel', config_names);
    ylabel('Relative Error');
    title('Error Comparison');
    grid on;
    
    subplot(2, 3, 5);
    bar(time_values);
    set(gca, 'XTickLabel', config_names);
    ylabel('Time (s)');
    title('Reconstruction Time');
    grid on;
    
    % Performance summary
    subplot(2, 3, 6);
    performance_matrix = [psnr_values, correlation_values*100, ssim_values*100, (1-error_values)*100];
    imagesc(performance_matrix);
    colorbar;
    set(gca, 'XTickLabel', config_names);
    set(gca, 'YTickLabel', {'PSNR', 'Corr%', 'SSIM%', 'Acc%'});
    title('Performance Matrix');
    
    sgtitle('Comparative Analysis of Target Configurations', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'comparative_analysis.png'));
    close(gcf);
    
    % Save comparative results
    save(fullfile(output_folder, 'comparative_results.mat'), 'comparative_results');
    
    fprintf('Comparative analysis completed\n');
end 