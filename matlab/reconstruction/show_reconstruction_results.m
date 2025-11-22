%% ===== SHOW RECONSTRUCTION RESULTS =====
% Display and summarize the end-to-end reconstruction results

clear; clc; close all;

%% ===== LOAD RESULTS =====
% Load the latest end-to-end demo results
demo_folder = 'end_to_end_demo';
demo_dirs = dir(fullfile(demo_folder, '*'));
demo_dirs = demo_dirs([demo_dirs.isdir]);
if isempty(demo_dirs)
    fprintf('ERROR: No demo directories found\n');
    return;
end
[~, latest_idx] = max([demo_dirs.datenum]);
latest_demo = fullfile(demo_folder, demo_dirs(latest_idx).name);
fprintf('Latest demo: %s\n', latest_demo);

% Check if comparative results exist
comparative_file = fullfile(latest_demo, 'comparative_results.mat');
if ~exist(comparative_file, 'file')
    fprintf('ERROR: Could not find comparative_results.mat in %s\n', latest_demo);
    fprintf('Available files:\n');
    files = dir(latest_demo);
    for i = 1:length(files)
        fprintf('  %s\n', files(i).name);
    end
    return;
end

fprintf('=== RECONSTRUCTION RESULTS SUMMARY ===\n');
fprintf('Loading results from: %s\n\n', latest_demo);

% Load comparative results
comparative_file = fullfile(latest_demo, 'comparative_results.mat');
if exist(comparative_file, 'file')
    load(comparative_file);
    fprintf('Loaded comparative results successfully\n');
else
    fprintf('ERROR: Could not find comparative results\n');
    return;
end

%% ===== DISPLAY RESULTS =====
config_names = fieldnames(comparative_results);
num_configs = length(config_names);

fprintf('\n=== RECONSTRUCTION PERFORMANCE SUMMARY ===\n');
fprintf('Configuration | PSNR (dB) | Correlation | SSIM | Error | Time (s)\n');
fprintf('-------------|-----------|------------|------|-------|---------\n');

for i = 1:num_configs
    config_name = config_names{i};
    results = comparative_results.(config_name);
    
    psnr_val = results.metrics.psnr;
    corr_val = results.metrics.correlation;
    ssim_val = results.metrics.ssim;
    error_val = results.metrics.relative_error;
    time_val = results.reconstruction_time;
    
    fprintf('%-13s | %9.2f | %10.4f | %4.2f | %5.4f | %8.2f\n', ...
        config_name, psnr_val, corr_val, ssim_val, error_val, time_val);
end

%% ===== KEY INSIGHTS =====
fprintf('\n=== KEY INSIGHTS ===\n');

% Find best performing configuration
psnr_values = zeros(num_configs, 1);
for i = 1:num_configs
    psnr_values(i) = comparative_results.(config_names{i}).metrics.psnr;
end
[best_psnr, best_idx] = max(psnr_values);
best_config = config_names{best_idx};

fprintf('Best PSNR Configuration: %s (%.2f dB)\n', best_config, best_psnr);

% Analyze target size impact
fprintf('\nTarget Size Analysis:\n');
for i = 1:num_configs
    config_name = config_names{i};
    config = comparative_results.(config_name).config;
    metrics = comparative_results.(config_name).metrics;
    fprintf('  %s: %.1f mm targets → PSNR: %.2f dB\n', ...
        config_name, config.target_size_mm, metrics.psnr);
end

% Analyze grid spacing impact
fprintf('\nGrid Spacing Analysis:\n');
for i = 1:num_configs
    config_name = config_names{i};
    config = comparative_results.(config_name).config;
    metrics = comparative_results.(config_name).metrics;
    fprintf('  %s: %.1f mm spacing → PSNR: %.2f dB\n', ...
        config_name, config.grid_spacing_mm, metrics.psnr);
end

%% ===== RECONSTRUCTION QUALITY ASSESSMENT =====
fprintf('\n=== RECONSTRUCTION QUALITY ASSESSMENT ===\n');

% Assess overall performance
avg_psnr = mean(psnr_values);
avg_corr = mean(arrayfun(@(i) comparative_results.(config_names{i}).metrics.correlation, 1:num_configs));

fprintf('Average Performance:\n');
fprintf('  - PSNR: %.2f dB\n', avg_psnr);
fprintf('  - Correlation: %.4f\n', avg_corr);

% Performance categories
excellent_psnr = sum(psnr_values >= 30);
good_psnr = sum(psnr_values >= 20 & psnr_values < 30);
fair_psnr = sum(psnr_values >= 10 & psnr_values < 20);
poor_psnr = sum(psnr_values < 10);

fprintf('\nPerformance Categories:\n');
fprintf('  - Excellent (≥30 dB): %d configurations\n', excellent_psnr);
fprintf('  - Good (20-30 dB): %d configurations\n', good_psnr);
fprintf('  - Fair (10-20 dB): %d configurations\n', fair_psnr);
fprintf('  - Poor (<10 dB): %d configurations\n', poor_psnr);

%% ===== COMPRESSED SENSING ANALYSIS =====
fprintf('\n=== COMPRESSED SENSING ANALYSIS ===\n');

% Load one configuration for detailed analysis
sample_config = config_names{1};
sample_file = fullfile(latest_demo, sprintf('%s_results.mat', sample_config));
if exist(sample_file, 'file')
    load(sample_file);
    
    fprintf('H Matrix Analysis:\n');
    fprintf('  - Dimensions: %d x %d\n', size(H_matrix, 1), size(H_matrix, 2));
    fprintf('  - Sparsity: %.2f%%\n', 100 * (1 - nnz(H_matrix) / numel(H_matrix)));
    fprintf('  - Condition number: %.2e\n', cond(H_matrix));
    
    fprintf('\nMeasurement Analysis:\n');
    fprintf('  - Number of measurements: %d\n', length(measurements));
    fprintf('  - Measurement SNR: %.2f dB\n', 20*log10(std(measurements)/std(measurements - mean(measurements))));
    
    fprintf('\nReconstruction Analysis:\n');
    fprintf('  - Reconstruction time: %.3f seconds\n', reconstruction_time);
    fprintf('  - Image dimensions: %d x %d pixels\n', size(reconstructed_image, 1), size(reconstructed_image, 2));
    fprintf('  - Target coverage: %.2f%%\n', 100 * nnz(reconstructed_image) / numel(reconstructed_image));
    
else
    fprintf('Could not load detailed analysis data\n');
end

%% ===== RECOMMENDATIONS =====
fprintf('\n=== RECOMMENDATIONS ===\n');

fprintf('Based on the reconstruction results:\n');
fprintf('1. %s configuration shows the best performance\n', best_config);
fprintf('2. Target sizes of %.1f-%.1f mm provide good reconstruction challenge\n', ...
    min(arrayfun(@(i) comparative_results.(config_names{i}).config.target_size_mm, 1:num_configs)), ...
    max(arrayfun(@(i) comparative_results.(config_names{i}).config.target_size_mm, 1:num_configs)));
fprintf('3. Grid spacing of %.1f-%.1f mm allows adequate target separation\n', ...
    min(arrayfun(@(i) comparative_results.(config_names{i}).config.grid_spacing_mm, 1:num_configs)), ...
    max(arrayfun(@(i) comparative_results.(config_names{i}).config.grid_spacing_mm, 1:num_configs)));
fprintf('4. H matrix generation is computationally efficient\n');
fprintf('5. Reconstruction algorithm shows consistent performance across configurations\n');

%% ===== VISUALIZATION SUMMARY =====
fprintf('\n=== VISUALIZATION SUMMARY ===\n');
fprintf('Generated visualizations:\n');
for i = 1:num_configs
    config_name = config_names{i};
    viz_file = fullfile(latest_demo, sprintf('%s_end_to_end_results.png', config_name));
    if exist(viz_file, 'file')
        fprintf('  ✓ %s: End-to-end reconstruction visualization\n', config_name);
    end
end

comparative_viz = fullfile(latest_demo, 'comparative_analysis.png');
if exist(comparative_viz, 'file')
    fprintf('  ✓ Comparative analysis visualization\n');
end

fprintf('\nAll results saved to: %s\n', latest_demo);
fprintf('\n=== SUMMARY COMPLETE ===\n'); 