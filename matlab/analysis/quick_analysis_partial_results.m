%% quick_analysis_partial_results.m - Quick Analysis of Partial Sweep Results
% Analyze the 17 completed runs from the interrupted parameter sweep

clear; clc; close all;

%% ===== CONFIGURATION =====
sweep_output_folder = 'realistic_parameter_sweep_output/072625_141142';

fprintf('=== QUICK ANALYSIS OF PARTIAL SWEEP RESULTS ===\n');
fprintf('Analyzing results from: %s\n\n', sweep_output_folder);

%% ===== EXTRACT RESULTS FROM COMPLETED RUNS =====
fprintf('Extracting results from completed runs...\n');

% Get all run folders
run_folders = dir(fullfile(sweep_output_folder, 'run*'));
num_runs = length(run_folders);

fprintf('Found %d completed runs\n', num_runs);

% Initialize arrays
psnr_values = zeros(num_runs, 1);
correlation_values = zeros(num_runs, 1);
reconstruction_times = zeros(num_runs, 1);
grid_steps = zeros(num_runs, 1);
target_sizes = zeros(num_runs, 1);
grid_spacings = zeros(num_runs, 1);
run_names = cell(num_runs, 1);

% Extract data from each run
for i = 1:num_runs
    run_folder = fullfile(sweep_output_folder, run_folders(i).name);
    results_file = fullfile(run_folder, 'run_results.mat');
    
    if exist(results_file, 'file')
        run_data = load(results_file);
        
        % Extract metrics
        psnr_values(i) = run_data.psnr;
        correlation_values(i) = run_data.correlation;
        reconstruction_times(i) = run_data.reconstruction_time;
        
        % Extract parameters from folder name
        folder_parts = strsplit(run_folders(i).name, '_');
        grid_steps(i) = str2double(folder_parts{2}(3:end));
        target_sizes(i) = str2double(folder_parts{3}(3:end));
        grid_spacings(i) = str2double(folder_parts{4}(3:end));
        run_names{i} = run_folders(i).name;
        
        fprintf('Run %d: PSNR=%.2f dB, Corr=%.4f, Time=%.1f s, GS=%.1f, TS=%.1f, SP=%.1f\n', ...
            i, psnr_values(i), correlation_values(i), reconstruction_times(i), ...
            grid_steps(i), target_sizes(i), grid_spacings(i));
    end
end

%% ===== BASIC STATISTICS =====
fprintf('\n=== BASIC STATISTICS (17 runs) ===\n');
fprintf('PSNR: %.2f ± %.2f dB (range: %.2f - %.2f dB)\n', ...
    mean(psnr_values), std(psnr_values), min(psnr_values), max(psnr_values));
fprintf('Correlation: %.4f ± %.4f (range: %.4f - %.4f)\n', ...
    mean(correlation_values), std(correlation_values), min(correlation_values), max(correlation_values));
fprintf('Average reconstruction time: %.1f ± %.1f seconds\n', ...
    mean(reconstruction_times), std(reconstruction_times));

%% ===== RANKING ANALYSIS =====
fprintf('\n=== RANKING BY PSNR ===\n');
[sorted_psnr, sort_indices] = sort(psnr_values, 'descend');

fprintf('Rank | PSNR (dB) | Correlation | Grid Step | Target Size | Grid Spacing | Time (s) | Run Name\n');
fprintf('-----|-----------|------------|----------|-------------|--------------|----------|---------\n');

for rank = 1:num_runs
    idx = sort_indices(rank);
    fprintf('%4d | %9.2f | %10.4f | %8.1f | %11.1f | %13.1f | %8.1f | %s\n', ...
        rank, psnr_values(idx), correlation_values(idx), grid_steps(idx), target_sizes(idx), ...
        grid_spacings(idx), reconstruction_times(idx), run_names{idx});
end

%% ===== PARAMETER ANALYSIS =====
fprintf('\n=== PARAMETER ANALYSIS ===\n');

% Grid step analysis
unique_gs = unique(grid_steps);
fprintf('Grid Step Analysis:\n');
for gs = unique_gs'
    mask = grid_steps == gs;
    avg_psnr = mean(psnr_values(mask));
    avg_corr = mean(correlation_values(mask));
    fprintf('  %.1f mm: PSNR=%.2f dB, Corr=%.4f (%d runs)\n', gs, avg_psnr, avg_corr, sum(mask));
end

% Target size analysis
unique_ts = unique(target_sizes);
fprintf('\nTarget Size Analysis:\n');
for ts = unique_ts'
    mask = target_sizes == ts;
    avg_psnr = mean(psnr_values(mask));
    avg_corr = mean(correlation_values(mask));
    fprintf('  %.1f mm: PSNR=%.2f dB, Corr=%.4f (%d runs)\n', ts, avg_psnr, avg_corr, sum(mask));
end

% Grid spacing analysis
unique_sp = unique(grid_spacings);
fprintf('\nGrid Spacing Analysis:\n');
for sp = unique_sp'
    mask = grid_spacings == sp;
    avg_psnr = mean(psnr_values(mask));
    avg_corr = mean(correlation_values(mask));
    fprintf('  %.1f mm: PSNR=%.2f dB, Corr=%.4f (%d runs)\n', sp, avg_psnr, avg_corr, sum(mask));
end

%% ===== BEST CONFIGURATION ===
best_idx = sort_indices(1);
fprintf('\n=== BEST CONFIGURATION ===\n');
fprintf('PSNR: %.2f dB\n', psnr_values(best_idx));
fprintf('Correlation: %.4f\n', correlation_values(best_idx));
fprintf('Grid Step: %.1f mm\n', grid_steps(best_idx));
fprintf('Target Size: %.1f mm\n', target_sizes(best_idx));
fprintf('Grid Spacing: %.1f mm\n', grid_spacings(best_idx));
fprintf('Reconstruction Time: %.1f seconds\n', reconstruction_times(best_idx));
fprintf('Run Name: %s\n', run_names{best_idx});

%% ===== QUICK VISUALIZATIONS ===
fprintf('\nCreating quick visualizations...\n');

% PSNR vs parameters
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
scatter(grid_steps, psnr_values, 100, psnr_values, 'filled');
colorbar;
xlabel('Grid Step (mm)');
ylabel('PSNR (dB)');
title('PSNR vs Grid Step');
grid on;

subplot(1, 3, 2);
scatter(target_sizes, psnr_values, 100, psnr_values, 'filled');
colorbar;
xlabel('Target Size (mm)');
ylabel('PSNR (dB)');
title('PSNR vs Target Size');
grid on;

subplot(1, 3, 3);
scatter(grid_spacings, psnr_values, 100, psnr_values, 'filled');
colorbar;
xlabel('Grid Spacing (mm)');
ylabel('PSNR (dB)');
title('PSNR vs Grid Spacing');
grid on;

sgtitle('Partial Parameter Sweep Results (17 runs)', 'FontSize', 14);

% Save plot
saveas(gcf, 'partial_sweep_analysis.png');
fprintf('Plot saved as: partial_sweep_analysis.png\n');

%% ===== KEY INSIGHTS ===
fprintf('\n=== KEY INSIGHTS ===\n');

% Best performing grid step
[~, best_gs_idx] = max(arrayfun(@(x) mean(psnr_values(grid_steps == x)), unique_gs));
fprintf('Best Grid Step: %.1f mm (avg PSNR: %.2f dB)\n', ...
    unique_gs(best_gs_idx), mean(psnr_values(grid_steps == unique_gs(best_gs_idx))));

% Best performing target size
[~, best_ts_idx] = max(arrayfun(@(x) mean(psnr_values(target_sizes == x)), unique_ts));
fprintf('Best Target Size: %.1f mm (avg PSNR: %.2f dB)\n', ...
    unique_ts(best_ts_idx), mean(psnr_values(target_sizes == unique_ts(best_ts_idx))));

% Best performing grid spacing
[~, best_sp_idx] = max(arrayfun(@(x) mean(psnr_values(grid_spacings == x)), unique_sp));
fprintf('Best Grid Spacing: %.1f mm (avg PSNR: %.2f dB)\n', ...
    unique_sp(best_sp_idx), mean(psnr_values(grid_spacings == unique_sp(best_sp_idx))));

% Performance categories
excellent_psnr = sum(psnr_values >= 40);
good_psnr = sum(psnr_values >= 30 & psnr_values < 40);
fair_psnr = sum(psnr_values >= 20 & psnr_values < 30);
poor_psnr = sum(psnr_values < 20);

fprintf('\nPerformance Categories:\n');
fprintf('  Excellent (≥40 dB): %d runs (%.1f%%)\n', excellent_psnr, excellent_psnr/num_runs*100);
fprintf('  Good (30-40 dB): %d runs (%.1f%%)\n', good_psnr, good_psnr/num_runs*100);
fprintf('  Fair (20-30 dB): %d runs (%.1f%%)\n', fair_psnr, fair_psnr/num_runs*100);
fprintf('  Poor (<20 dB): %d runs (%.1f%%)\n', poor_psnr, poor_psnr/num_runs*100);

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Best PSNR achieved: %.2f dB\n', max(psnr_values));
fprintf('Best correlation achieved: %.4f\n', max(correlation_values)); 