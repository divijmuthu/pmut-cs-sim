%% analyze_parameter_sweep_results.m - Comprehensive Parameter Sweep Analysis
% Analyzes results from realistic_parameter_sweep.m and generates insights
% Run this after the parameter sweep completes

clear; clc; close all;

%% ===== CONFIGURATION =====
% Point to the sweep output folder (update timestamp as needed)
sweep_output_folder = 'realistic_parameter_sweep_output/072625_141142'; % Update this timestamp

% Analysis output folder
analysis_timestamp = datestr(now, 'mmddyy_HHMMSS');
analysis_output_folder = fullfile('sweep_analysis_output', analysis_timestamp);
if ~exist(analysis_output_folder, 'dir')
    mkdir(analysis_output_folder);
end

fprintf('=== PARAMETER SWEEP ANALYSIS ===\n');
fprintf('Analyzing results from: %s\n', sweep_output_folder);
fprintf('Saving analysis to: %s\n\n', analysis_output_folder);

%% ===== LOAD SWEEP RESULTS =====
fprintf('Loading sweep results...\n');

% Load the complete results
results_file = fullfile(sweep_output_folder, 'complete_sweep_results.mat');
if exist(results_file, 'file')
    load(results_file);
    fprintf('Loaded complete results from: %s\n', results_file);
else
    fprintf('ERROR: Could not find complete_sweep_results.mat\n');
    fprintf('Checking for individual run folders...\n');
    
    % Try to reconstruct results from individual run folders
    [results, param_combinations] = reconstruct_results_from_folders(sweep_output_folder);
    if isempty(results)
        error('Could not load or reconstruct sweep results');
    end
end

num_runs = length(results.metrics);
fprintf('Found %d parameter combinations\n\n', num_runs);

%% ===== EXTRACT METRICS =====
fprintf('Extracting metrics...\n');

% Extract all metrics into arrays
psnr_values = zeros(num_runs, 1);
correlation_values = zeros(num_runs, 1);
reconstruction_times = zeros(num_runs, 1);

% Extract parameters
grid_steps = zeros(num_runs, 1);
target_sizes = zeros(num_runs, 1);
grid_spacings = zeros(num_runs, 1);
num_acquisitions = zeros(num_runs, 1);
lambda_tv_values = zeros(num_runs, 1);
target_patterns = cell(num_runs, 1);

for i = 1:num_runs
    if ~isempty(results.metrics{i})
        psnr_values(i) = results.metrics{i}.psnr;
        correlation_values(i) = results.metrics{i}.correlation;
        reconstruction_times(i) = results.times{i}.total_time;
        
        config = results.configs{i};
        grid_steps(i) = config.grid_step_mm;
        target_sizes(i) = config.target_size_mm;
        grid_spacings(i) = config.grid_spacing_mm;
        num_acquisitions(i) = config.num_acquisitions;
        lambda_tv_values(i) = config.lambda_tv_reg;
        target_patterns{i} = config.target_pattern;
    end
end

%% ===== BASIC STATISTICS =====
fprintf('Computing basic statistics...\n');

% PSNR statistics
psnr_stats = struct();
psnr_stats.mean = mean(psnr_values);
psnr_stats.std = std(psnr_values);
psnr_stats.min = min(psnr_values);
psnr_stats.max = max(psnr_values);
psnr_stats.median = median(psnr_values);

% Correlation statistics
corr_stats = struct();
corr_stats.mean = mean(correlation_values);
corr_stats.std = std(correlation_values);
corr_stats.min = min(correlation_values);
corr_stats.max = max(correlation_values);

% Time statistics
time_stats = struct();
time_stats.mean = mean(reconstruction_times);
time_stats.std = std(reconstruction_times);
time_stats.total = sum(reconstruction_times);

fprintf('\n=== BASIC STATISTICS ===\n');
fprintf('PSNR: %.2f ± %.2f dB (range: %.2f - %.2f dB)\n', ...
    psnr_stats.mean, psnr_stats.std, psnr_stats.min, psnr_stats.max);
fprintf('Correlation: %.4f ± %.4f (range: %.4f - %.4f)\n', ...
    corr_stats.mean, corr_stats.std, corr_stats.min, corr_stats.max);
fprintf('Average reconstruction time: %.1f ± %.1f seconds\n', ...
    time_stats.mean, time_stats.std);
fprintf('Total sweep time: %.1f minutes\n', time_stats.total / 60);

%% ===== RANKING ANALYSIS =====
fprintf('\nRanking results by PSNR...\n');

% Sort by PSNR (descending)
[sorted_psnr, sort_indices] = sort(psnr_values, 'descend');

% Top 10 results
fprintf('\n=== TOP 10 RESULTS BY PSNR ===\n');
fprintf('Rank | PSNR (dB) | Correlation | Grid Step | Target Size | Grid Spacing | Acquisitions | Lambda TV | Pattern | Time (s)\n');
fprintf('-----|-----------|------------|----------|-------------|--------------|-------------|-----------|---------|---------\n');

for rank = 1:min(10, num_runs)
    idx = sort_indices(rank);
    fprintf('%4d | %9.2f | %10.4f | %8.1f | %11.1f | %13.1f | %12d | %9.4f | %7s | %8.1f\n', ...
        rank, psnr_values(idx), correlation_values(idx), grid_steps(idx), target_sizes(idx), ...
        grid_spacings(idx), num_acquisitions(idx), lambda_tv_values(idx), ...
        target_patterns{idx}, reconstruction_times(idx));
end

% Bottom 10 results
fprintf('\n=== BOTTOM 10 RESULTS BY PSNR ===\n');
fprintf('Rank | PSNR (dB) | Correlation | Grid Step | Target Size | Grid Spacing | Acquisitions | Lambda TV | Pattern | Time (s)\n');
fprintf('-----|-----------|------------|----------|-------------|--------------|-------------|-----------|---------|---------\n');

for rank = max(1, num_runs-9):num_runs
    idx = sort_indices(rank);
    fprintf('%4d | %9.2f | %10.4f | %8.1f | %11.1f | %13.1f | %12d | %9.4f | %7s | %8.1f\n', ...
        rank, psnr_values(idx), correlation_values(idx), grid_steps(idx), target_sizes(idx), ...
        grid_spacings(idx), num_acquisitions(idx), lambda_tv_values(idx), ...
        target_patterns{idx}, reconstruction_times(idx));
end

%% ===== PARAMETER SENSITIVITY ANALYSIS =====
fprintf('\nAnalyzing parameter sensitivity...\n');

% Create parameter sensitivity plots
create_parameter_sensitivity_plots(psnr_values, correlation_values, reconstruction_times, ...
    grid_steps, target_sizes, grid_spacings, num_acquisitions, lambda_tv_values, ...
    target_patterns, analysis_output_folder);

%% ===== BEST CONFIGURATION ANALYSIS =====
fprintf('\nAnalyzing best configuration...\n');

% Get best configuration
best_idx = sort_indices(1);
best_config = results.configs{best_idx};
best_psnr = psnr_values(best_idx);
best_correlation = correlation_values(best_idx);

fprintf('\n=== BEST CONFIGURATION ===\n');
fprintf('PSNR: %.2f dB\n', best_psnr);
fprintf('Correlation: %.4f\n', best_correlation);
fprintf('Grid Step: %.1f mm\n', best_config.grid_step_mm);
fprintf('Target Size: %.1f mm\n', best_config.target_size_mm);
fprintf('Grid Spacing: %.1f mm\n', best_config.grid_spacing_mm);
fprintf('Acquisitions: %d\n', best_config.num_acquisitions);
fprintf('Lambda TV: %.4f\n', best_config.lambda_tv_reg);
fprintf('Pattern: %s\n', best_config.target_pattern);

% Save best configuration
save(fullfile(analysis_output_folder, 'best_configuration.mat'), 'best_config', 'best_psnr', 'best_correlation');

%% ===== PERFORMANCE CATEGORIES =====
fprintf('\nCategorizing performance...\n');

% Define performance categories
excellent_psnr = psnr_values >= 40;  % Excellent: PSNR >= 40 dB
good_psnr = psnr_values >= 30 & psnr_values < 40;  % Good: 30-40 dB
fair_psnr = psnr_values >= 20 & psnr_values < 30;  % Fair: 20-30 dB
poor_psnr = psnr_values < 20;  % Poor: < 20 dB

excellent_corr = correlation_values >= 0.99;  % Excellent: Corr >= 0.99
good_corr = correlation_values >= 0.95 & correlation_values < 0.99;  % Good: 0.95-0.99
fair_corr = correlation_values >= 0.90 & correlation_values < 0.95;  % Fair: 0.90-0.95
poor_corr = correlation_values < 0.90;  % Poor: < 0.90

fprintf('\n=== PERFORMANCE CATEGORIES ===\n');
fprintf('PSNR Categories:\n');
fprintf('  Excellent (≥40 dB): %d runs (%.1f%%)\n', sum(excellent_psnr), sum(excellent_psnr)/num_runs*100);
fprintf('  Good (30-40 dB): %d runs (%.1f%%)\n', sum(good_psnr), sum(good_psnr)/num_runs*100);
fprintf('  Fair (20-30 dB): %d runs (%.1f%%)\n', sum(fair_psnr), sum(fair_psnr)/num_runs*100);
fprintf('  Poor (<20 dB): %d runs (%.1f%%)\n', sum(poor_psnr), sum(poor_psnr)/num_runs*100);

fprintf('\nCorrelation Categories:\n');
fprintf('  Excellent (≥0.99): %d runs (%.1f%%)\n', sum(excellent_corr), sum(excellent_corr)/num_runs*100);
fprintf('  Good (0.95-0.99): %d runs (%.1f%%)\n', sum(good_corr), sum(good_corr)/num_runs*100);
fprintf('  Fair (0.90-0.95): %d runs (%.1f%%)\n', sum(fair_corr), sum(fair_corr)/num_runs*100);
fprintf('  Poor (<0.90): %d runs (%.1f%%)\n', sum(poor_corr), sum(poor_corr)/num_runs*100);

%% ===== PARAMETER OPTIMIZATION INSIGHTS =====
fprintf('\nGenerating optimization insights...\n');

% Analyze parameter trends
create_optimization_insights(psnr_values, correlation_values, ...
    grid_steps, target_sizes, grid_spacings, num_acquisitions, lambda_tv_values, ...
    target_patterns, analysis_output_folder);

%% ===== COMPREHENSIVE VISUALIZATION =====
fprintf('\nCreating comprehensive visualizations...\n');

% Create summary dashboard
create_summary_dashboard(psnr_values, correlation_values, reconstruction_times, ...
    grid_steps, target_sizes, grid_spacings, num_acquisitions, lambda_tv_values, ...
    target_patterns, sort_indices, analysis_output_folder);

%% ===== SAVE ANALYSIS RESULTS =====
fprintf('\nSaving analysis results...\n');

% Save comprehensive analysis
analysis_results = struct();
analysis_results.psnr_stats = psnr_stats;
analysis_results.corr_stats = corr_stats;
analysis_results.time_stats = time_stats;
analysis_results.best_config = best_config;
analysis_results.best_psnr = best_psnr;
analysis_results.best_correlation = best_correlation;
analysis_results.sort_indices = sort_indices;
analysis_results.performance_categories = struct();
analysis_results.performance_categories.excellent_psnr = sum(excellent_psnr);
analysis_results.performance_categories.good_psnr = sum(good_psnr);
analysis_results.performance_categories.fair_psnr = sum(fair_psnr);
analysis_results.performance_categories.poor_psnr = sum(poor_psnr);

save(fullfile(analysis_output_folder, 'comprehensive_analysis.mat'), 'analysis_results');

% Create detailed CSV report
create_detailed_csv_report(results, sort_indices, analysis_output_folder);

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Results saved to: %s\n', analysis_output_folder);
fprintf('Best PSNR achieved: %.2f dB\n', best_psnr);
fprintf('Best correlation achieved: %.4f\n', best_correlation);

%% ===== HELPER FUNCTIONS =====

function [results, param_combinations] = reconstruct_results_from_folders(sweep_folder)
    % Reconstruct results from individual run folders if complete results not available
    
    fprintf('Reconstructing results from individual run folders...\n');
    
    % Get all run folders
    run_folders = dir(fullfile(sweep_folder, 'run*'));
    num_runs = length(run_folders);
    
    if num_runs == 0
        results = [];
        param_combinations = [];
        return;
    end
    
    % Initialize results structure
    results = struct();
    results.configs = cell(num_runs, 1);
    results.metrics = cell(num_runs, 1);
    results.times = cell(num_runs, 1);
    results.psnr_rankings = zeros(num_runs, 1);
    
    param_combinations = [];
    
    for i = 1:num_runs
        run_folder = fullfile(sweep_folder, run_folders(i).name);
        results_file = fullfile(run_folder, 'run_results.mat');
        
        if exist(results_file, 'file')
            run_data = load(results_file);
            
            % Extract parameters from folder name
            folder_parts = strsplit(run_folders(i).name, '_');
            
            config = struct();
            config.grid_step_mm = str2double(folder_parts{2}(3:end));
            config.target_size_mm = str2double(folder_parts{3}(3:end));
            config.grid_spacing_mm = str2double(folder_parts{4}(3:end));
            config.num_acquisitions = str2double(folder_parts{5}(4:end));
            config.lambda_tv_reg = str2double(folder_parts{6}(4:end));
            config.target_pattern = folder_parts{7};
            
            metrics = struct();
            metrics.psnr = run_data.psnr;
            metrics.correlation = run_data.correlation;
            
            times = struct();
            times.total_time = run_data.reconstruction_time;
            
            results.configs{i} = config;
            results.metrics{i} = metrics;
            results.times{i} = times;
            results.psnr_rankings(i) = run_data.psnr;
            
            param_combinations = [param_combinations; config];
        end
    end
    
    fprintf('Reconstructed results for %d runs\n', num_runs);
end

function create_parameter_sensitivity_plots(psnr_values, correlation_values, reconstruction_times, ...
    grid_steps, target_sizes, grid_spacings, num_acquisitions, lambda_tv_values, ...
    target_patterns, output_folder)
    
    % Create parameter sensitivity plots
    figure('Position', [100, 100, 1200, 800]);
    
    % PSNR vs parameters
    subplot(2, 3, 1);
    scatter(grid_steps, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Step (mm)');
    ylabel('PSNR (dB)');
    title('PSNR vs Grid Step');
    grid on;
    
    subplot(2, 3, 2);
    scatter(target_sizes, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Target Size (mm)');
    ylabel('PSNR (dB)');
    title('PSNR vs Target Size');
    grid on;
    
    subplot(2, 3, 3);
    scatter(grid_spacings, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Spacing (mm)');
    ylabel('PSNR (dB)');
    title('PSNR vs Grid Spacing');
    grid on;
    
    subplot(2, 3, 4);
    scatter(num_acquisitions, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Number of Acquisitions');
    ylabel('PSNR (dB)');
    title('PSNR vs Acquisitions');
    grid on;
    
    subplot(2, 3, 5);
    scatter(lambda_tv_values, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Lambda TV');
    ylabel('PSNR (dB)');
    title('PSNR vs Lambda TV');
    grid on;
    
    subplot(2, 3, 6);
    % Pattern analysis
    unique_patterns = unique(target_patterns);
    pattern_psnr = zeros(length(unique_patterns), 1);
    for i = 1:length(unique_patterns)
        pattern_mask = strcmp(target_patterns, unique_patterns{i});
        pattern_psnr(i) = mean(psnr_values(pattern_mask));
    end
    bar(pattern_psnr);
    set(gca, 'XTickLabel', unique_patterns);
    xlabel('Target Pattern');
    ylabel('Average PSNR (dB)');
    title('PSNR vs Target Pattern');
    grid on;
    
    saveas(gcf, fullfile(output_folder, 'parameter_sensitivity_psnr.png'));
    close(gcf);
    
    % Correlation vs parameters
    figure('Position', [100, 100, 1200, 800]);
    
    subplot(2, 3, 1);
    scatter(grid_steps, correlation_values, 50, correlation_values, 'filled');
    colorbar;
    xlabel('Grid Step (mm)');
    ylabel('Correlation');
    title('Correlation vs Grid Step');
    grid on;
    
    subplot(2, 3, 2);
    scatter(target_sizes, correlation_values, 50, correlation_values, 'filled');
    colorbar;
    xlabel('Target Size (mm)');
    ylabel('Correlation');
    title('Correlation vs Target Size');
    grid on;
    
    subplot(2, 3, 3);
    scatter(grid_spacings, correlation_values, 50, correlation_values, 'filled');
    colorbar;
    xlabel('Grid Spacing (mm)');
    ylabel('Correlation');
    title('Correlation vs Grid Spacing');
    grid on;
    
    subplot(2, 3, 4);
    scatter(num_acquisitions, correlation_values, 50, correlation_values, 'filled');
    colorbar;
    xlabel('Number of Acquisitions');
    ylabel('Correlation');
    title('Correlation vs Acquisitions');
    grid on;
    
    subplot(2, 3, 5);
    scatter(lambda_tv_values, correlation_values, 50, correlation_values, 'filled');
    colorbar;
    xlabel('Lambda TV');
    ylabel('Correlation');
    title('Correlation vs Lambda TV');
    grid on;
    
    subplot(2, 3, 6);
    % Pattern analysis
    pattern_corr = zeros(length(unique_patterns), 1);
    for i = 1:length(unique_patterns)
        pattern_mask = strcmp(target_patterns, unique_patterns{i});
        pattern_corr(i) = mean(correlation_values(pattern_mask));
    end
    bar(pattern_corr);
    set(gca, 'XTickLabel', unique_patterns);
    xlabel('Target Pattern');
    ylabel('Average Correlation');
    title('Correlation vs Target Pattern');
    grid on;
    
    saveas(gcf, fullfile(output_folder, 'parameter_sensitivity_correlation.png'));
    close(gcf);
end

function create_optimization_insights(psnr_values, correlation_values, ...
    grid_steps, target_sizes, grid_spacings, num_acquisitions, lambda_tv_values, ...
    target_patterns, output_folder)
    
    % Create optimization insights report
    insights_file = fullfile(output_folder, 'optimization_insights.txt');
    fid = fopen(insights_file, 'w');
    
    fprintf(fid, '=== PARAMETER SWEEP OPTIMIZATION INSIGHTS ===\n\n');
    
    % Grid step insights
    unique_grid_steps = unique(grid_steps);
    fprintf(fid, 'GRID STEP ANALYSIS:\n');
    for gs = unique_grid_steps'
        mask = grid_steps == gs;
        avg_psnr = mean(psnr_values(mask));
        avg_corr = mean(correlation_values(mask));
        fprintf(fid, '  %.1f mm: PSNR=%.2f dB, Corr=%.4f\n', gs, avg_psnr, avg_corr);
    end
    [~, best_gs_idx] = max(arrayfun(@(x) mean(psnr_values(grid_steps == x)), unique_grid_steps));
    fprintf(fid, '  RECOMMENDATION: %.1f mm grid step\n\n', unique_grid_steps(best_gs_idx));
    
    % Target size insights
    unique_target_sizes = unique(target_sizes);
    fprintf(fid, 'TARGET SIZE ANALYSIS:\n');
    for ts = unique_target_sizes'
        mask = target_sizes == ts;
        avg_psnr = mean(psnr_values(mask));
        avg_corr = mean(correlation_values(mask));
        fprintf(fid, '  %.1f mm: PSNR=%.2f dB, Corr=%.4f\n', ts, avg_psnr, avg_corr);
    end
    [~, best_ts_idx] = max(arrayfun(@(x) mean(psnr_values(target_sizes == x)), unique_target_sizes));
    fprintf(fid, '  RECOMMENDATION: %.1f mm target size\n\n', unique_target_sizes(best_ts_idx));
    
    % Grid spacing insights
    unique_grid_spacings = unique(grid_spacings);
    fprintf(fid, 'GRID SPACING ANALYSIS:\n');
    for gsp = unique_grid_spacings'
        mask = grid_spacings == gsp;
        avg_psnr = mean(psnr_values(mask));
        avg_corr = mean(correlation_values(mask));
        fprintf(fid, '  %.1f mm: PSNR=%.2f dB, Corr=%.4f\n', gsp, avg_psnr, avg_corr);
    end
    [~, best_gsp_idx] = max(arrayfun(@(x) mean(psnr_values(grid_spacings == x)), unique_grid_spacings));
    fprintf(fid, '  RECOMMENDATION: %.1f mm grid spacing\n\n', unique_grid_spacings(best_gsp_idx));
    
    % Acquisitions insights
    unique_acq = unique(num_acquisitions);
    fprintf(fid, 'ACQUISITIONS ANALYSIS:\n');
    for acq = unique_acq'
        mask = num_acquisitions == acq;
        avg_psnr = mean(psnr_values(mask));
        avg_corr = mean(correlation_values(mask));
        fprintf(fid, '  %d: PSNR=%.2f dB, Corr=%.4f\n', acq, avg_psnr, avg_corr);
    end
    [~, best_acq_idx] = max(arrayfun(@(x) mean(psnr_values(num_acquisitions == x)), unique_acq));
    fprintf(fid, '  RECOMMENDATION: %d acquisitions\n\n', unique_acq(best_acq_idx));
    
    % Lambda TV insights
    unique_lambda = unique(lambda_tv_values);
    fprintf(fid, 'LAMBDA TV ANALYSIS:\n');
    for lam = unique_lambda'
        mask = lambda_tv_values == lam;
        avg_psnr = mean(psnr_values(mask));
        avg_corr = mean(correlation_values(mask));
        fprintf(fid, '  %.4f: PSNR=%.2f dB, Corr=%.4f\n', lam, avg_psnr, avg_corr);
    end
    [~, best_lam_idx] = max(arrayfun(@(x) mean(psnr_values(lambda_tv_values == x)), unique_lambda));
    fprintf(fid, '  RECOMMENDATION: %.4f lambda TV\n\n', unique_lambda(best_lam_idx));
    
    % Pattern insights
    unique_patterns = unique(target_patterns);
    fprintf(fid, 'TARGET PATTERN ANALYSIS:\n');
    for i = 1:length(unique_patterns)
        pattern = unique_patterns{i};
        mask = strcmp(target_patterns, pattern);
        avg_psnr = mean(psnr_values(mask));
        avg_corr = mean(correlation_values(mask));
        fprintf(fid, '  %s: PSNR=%.2f dB, Corr=%.4f\n', pattern, avg_psnr, avg_corr);
    end
    [~, best_pat_idx] = max(arrayfun(@(x) mean(psnr_values(strcmp(target_patterns, x))), unique_patterns));
    fprintf(fid, '  RECOMMENDATION: %s pattern\n\n', unique_patterns{best_pat_idx});
    
    % Overall recommendations
    fprintf(fid, '=== OVERALL RECOMMENDATIONS ===\n');
    fprintf(fid, 'Best Grid Step: %.1f mm\n', unique_grid_steps(best_gs_idx));
    fprintf(fid, 'Best Target Size: %.1f mm\n', unique_target_sizes(best_ts_idx));
    fprintf(fid, 'Best Grid Spacing: %.1f mm\n', unique_grid_spacings(best_gsp_idx));
    fprintf(fid, 'Best Acquisitions: %d\n', unique_acq(best_acq_idx));
    fprintf(fid, 'Best Lambda TV: %.4f\n', unique_lambda(best_lam_idx));
    fprintf(fid, 'Best Pattern: %s\n', unique_patterns{best_pat_idx});
    
    fclose(fid);
    fprintf('Optimization insights saved to: %s\n', insights_file);
end

function create_summary_dashboard(psnr_values, correlation_values, reconstruction_times, ...
    grid_steps, target_sizes, grid_spacings, num_acquisitions, lambda_tv_values, ...
    target_patterns, sort_indices, output_folder)
    
    % Create comprehensive summary dashboard
    figure('Position', [50, 50, 1400, 1000]);
    
    % PSNR distribution
    subplot(3, 4, 1);
    histogram(psnr_values, 20, 'Normalization', 'probability');
    xlabel('PSNR (dB)');
    ylabel('Probability');
    title('PSNR Distribution');
    grid on;
    
    % Correlation distribution
    subplot(3, 4, 2);
    histogram(correlation_values, 20, 'Normalization', 'probability');
    xlabel('Correlation');
    ylabel('Probability');
    title('Correlation Distribution');
    grid on;
    
    % PSNR vs Correlation
    subplot(3, 4, 3);
    scatter(correlation_values, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Correlation');
    ylabel('PSNR (dB)');
    title('PSNR vs Correlation');
    grid on;
    
    % Time vs PSNR
    subplot(3, 4, 4);
    scatter(reconstruction_times, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Reconstruction Time (s)');
    ylabel('PSNR (dB)');
    title('PSNR vs Time');
    grid on;
    
    % Top 10 PSNR values
    subplot(3, 4, 5);
    top_10_psnr = psnr_values(sort_indices(1:min(10, length(sort_indices))));
    bar(top_10_psnr);
    xlabel('Rank');
    ylabel('PSNR (dB)');
    title('Top 10 PSNR Values');
    grid on;
    
    % Parameter heatmaps
    subplot(3, 4, 6);
    % Grid step vs Target size heatmap
    unique_gs = unique(grid_steps);
    unique_ts = unique(target_sizes);
    heatmap_data = zeros(length(unique_ts), length(unique_gs));
    for i = 1:length(unique_ts)
        for j = 1:length(unique_gs)
            mask = grid_steps == unique_gs(j) & target_sizes == unique_ts(i);
            if any(mask)
                heatmap_data(i, j) = mean(psnr_values(mask));
            end
        end
    end
    imagesc(heatmap_data);
    colorbar;
    set(gca, 'XTick', 1:length(unique_gs), 'XTickLabel', unique_gs);
    set(gca, 'YTick', 1:length(unique_ts), 'YTickLabel', unique_ts);
    xlabel('Grid Step (mm)');
    ylabel('Target Size (mm)');
    title('PSNR: Grid Step vs Target Size');
    
    subplot(3, 4, 7);
    % Grid spacing vs Acquisitions heatmap
    unique_gsp = unique(grid_spacings);
    unique_acq = unique(num_acquisitions);
    heatmap_data2 = zeros(length(unique_acq), length(unique_gsp));
    for i = 1:length(unique_acq)
        for j = 1:length(unique_gsp)
            mask = grid_spacings == unique_gsp(j) & num_acquisitions == unique_acq(i);
            if any(mask)
                heatmap_data2(i, j) = mean(psnr_values(mask));
            end
        end
    end
    imagesc(heatmap_data2);
    colorbar;
    set(gca, 'XTick', 1:length(unique_gsp), 'XTickLabel', unique_gsp);
    set(gca, 'YTick', 1:length(unique_acq), 'YTickLabel', unique_acq);
    xlabel('Grid Spacing (mm)');
    ylabel('Acquisitions');
    title('PSNR: Grid Spacing vs Acquisitions');
    
    subplot(3, 4, 8);
    % Lambda TV vs Pattern
    unique_lambda = unique(lambda_tv_values);
    unique_patterns = unique(target_patterns);
    heatmap_data3 = zeros(length(unique_patterns), length(unique_lambda));
    for i = 1:length(unique_patterns)
        for j = 1:length(unique_lambda)
            mask = strcmp(target_patterns, unique_patterns{i}) & lambda_tv_values == unique_lambda(j);
            if any(mask)
                heatmap_data3(i, j) = mean(psnr_values(mask));
            end
        end
    end
    imagesc(heatmap_data3);
    colorbar;
    set(gca, 'XTick', 1:length(unique_lambda), 'XTickLabel', unique_lambda);
    set(gca, 'YTick', 1:length(unique_patterns), 'YTickLabel', unique_patterns);
    xlabel('Lambda TV');
    ylabel('Pattern');
    title('PSNR: Lambda TV vs Pattern');
    
    % Performance categories
    subplot(3, 4, 9);
    excellent_psnr = sum(psnr_values >= 40);
    good_psnr = sum(psnr_values >= 30 & psnr_values < 40);
    fair_psnr = sum(psnr_values >= 20 & psnr_values < 30);
    poor_psnr = sum(psnr_values < 20);
    pie([excellent_psnr, good_psnr, fair_psnr, poor_psnr], ...
        {'Excellent (≥40 dB)', 'Good (30-40 dB)', 'Fair (20-30 dB)', 'Poor (<20 dB)'});
    title('PSNR Performance Categories');
    
    subplot(3, 4, 10);
    excellent_corr = sum(correlation_values >= 0.99);
    good_corr = sum(correlation_values >= 0.95 & correlation_values < 0.99);
    fair_corr = sum(correlation_values >= 0.90 & correlation_values < 0.95);
    poor_corr = sum(correlation_values < 0.90);
    pie([excellent_corr, good_corr, fair_corr, poor_corr], ...
        {'Excellent (≥0.99)', 'Good (0.95-0.99)', 'Fair (0.90-0.95)', 'Poor (<0.90)'});
    title('Correlation Performance Categories');
    
    % Parameter ranges
    subplot(3, 4, 11);
    param_ranges = [range(grid_steps), range(target_sizes), range(grid_spacings), ...
                   range(num_acquisitions), range(lambda_tv_values)];
    param_names = {'Grid Step', 'Target Size', 'Grid Spacing', 'Acquisitions', 'Lambda TV'};
    bar(param_ranges);
    set(gca, 'XTickLabel', param_names);
    ylabel('Range');
    title('Parameter Ranges Tested');
    grid on;
    
    % Best configuration summary
    subplot(3, 4, 12);
    best_idx = sort_indices(1);
    best_params = [grid_steps(best_idx), target_sizes(best_idx), grid_spacings(best_idx), ...
                   num_acquisitions(best_idx), lambda_tv_values(best_idx)];
    bar(best_params);
    set(gca, 'XTickLabel', param_names);
    ylabel('Value');
    title('Best Configuration Parameters');
    grid on;
    
    sgtitle('Parameter Sweep Analysis Dashboard', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'summary_dashboard.png'));
    close(gcf);
end

function create_detailed_csv_report(results, sort_indices, output_folder)
    % Create detailed CSV report with all results
    
    num_runs = length(results.metrics);
    report_data = [];
    
    for rank = 1:num_runs
        idx = sort_indices(rank);
        config = results.configs{idx};
        metrics = results.metrics{idx};
        times = results.times{idx};
        
        row = struct();
        row.rank = rank;
        row.psnr = metrics.psnr;
        row.correlation = metrics.correlation;
        row.reconstruction_time = times.total_time;
        row.grid_step_mm = config.grid_step_mm;
        row.target_size_mm = config.target_size_mm;
        row.grid_spacing_mm = config.grid_spacing_mm;
        row.num_acquisitions = config.num_acquisitions;
        row.lambda_tv_reg = config.lambda_tv_reg;
        row.target_pattern = config.target_pattern;
        
        report_data = [report_data; row];
    end
    
    % Convert to table and save
    report_table = struct2table(report_data);
    writetable(report_table, fullfile(output_folder, 'detailed_results_ranked.csv'));
    
    fprintf('Detailed CSV report saved to: %s\n', fullfile(output_folder, 'detailed_results_ranked.csv'));
end 