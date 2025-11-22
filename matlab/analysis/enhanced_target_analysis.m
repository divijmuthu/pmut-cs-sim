%% ===== ENHANCED TARGET ANALYSIS & VISUALIZATION =====
% Improved target setup based on TheoryTest2 with comprehensive visualization
% Designed for challenging reconstruction with 9 small, clearly distinguishable targets

clear; clc; close all;

%% ===== CONFIGURATION =====
% Analysis parameters
sweep_output_folder = 'realistic_parameter_sweep_output/072625_141142'; % Update timestamp
analysis_timestamp = datestr(now, 'mmddyy_HHMMSS');
analysis_output_folder = fullfile('enhanced_target_analysis', analysis_timestamp);
if ~exist(analysis_output_folder, 'dir')
    mkdir(analysis_output_folder);
end

fprintf('=== ENHANCED TARGET ANALYSIS ===\n');
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

%% ===== ENHANCED TARGET VISUALIZATION =====
fprintf('Creating enhanced target visualizations...\n');

% Create improved target scene based on TheoryTest2 approach
create_enhanced_target_scene(analysis_output_folder);

% Visualize sample reconstructions
visualize_sample_reconstructions(results, analysis_output_folder);

% Create target quality assessment
assess_target_quality(results, analysis_output_folder);

%% ===== IMPROVED TARGET SETUP =====
fprintf('Creating improved target setup recommendations...\n');

% Generate improved target configurations
improved_targets = create_improved_target_configurations();

% Save improved target configurations
save(fullfile(analysis_output_folder, 'improved_target_configurations.mat'), 'improved_targets');

%% ===== COMPREHENSIVE VISUALIZATION DASHBOARD =====
fprintf('Creating comprehensive visualization dashboard...\n');

create_comprehensive_dashboard(results, analysis_output_folder);

fprintf('\n=== ENHANCED ANALYSIS COMPLETE ===\n');
fprintf('Results saved to: %s\n', analysis_output_folder);

%% ===== HELPER FUNCTIONS =====

function create_enhanced_target_scene(output_folder)
    % Create enhanced target scene based on TheoryTest2 approach
    
    % Grid parameters (based on TheoryTest2)
    grid_width_mm = 150;
    grid_step_mm = 2;  % Higher resolution for better detail
    target_distance_mm = 250;
    
    % Create imaging grid
    x_coords = -grid_width_mm/2 : grid_step_mm : grid_width_mm/2;
    y_coords = -grid_width_mm/2 : grid_step_mm : grid_width_mm/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords, y_coords);
    
    % Enhanced target positions (based on TheoryTest2's 3x3 grid)
    target_positions = [
        -30, 30,  1.0;   % Top-left
        0,   30,  1.0;   % Top-center
        30,  30,  1.0;   % Top-right
        -30, 0,   1.0;   % Middle-left
        0,   0,   1.0;   % Center
        30,  0,   1.0;   % Middle-right
        -30, -30, 1.0;   % Bottom-left
        0,   -30, 1.0;   % Bottom-center
        30,  -30, 1.0    % Bottom-right
    ];
    
    % Create scene matrix
    scene_matrix = zeros(length(y_coords), length(x_coords));
    
    % Enhanced target creation with different sizes and intensities
    target_sizes_mm = [3, 4, 5, 6, 8];  % Different target sizes for challenge
    target_intensities = [0.7, 0.8, 0.9, 1.0, 1.1];  % Different intensities
    
    for i = 1:size(target_positions, 1)
        x_pos_mm = target_positions(i, 1);
        y_pos_mm = target_positions(i, 2);
        base_amplitude = target_positions(i, 3);
        
        % Vary target properties for more challenge
        target_size_mm = target_sizes_mm(mod(i-1, length(target_sizes_mm)) + 1);
        target_intensity = target_intensities(mod(i-1, length(target_intensities)) + 1);
        
        % Find nearest grid point
        [~, ix_scene] = min(abs(x_coords - x_pos_mm));
        [~, iy_scene] = min(abs(y_coords - y_pos_mm));
        
        % Create target with specified size
        target_radius_pixels = round(target_size_mm / (2 * grid_step_mm));
        if target_radius_pixels < 1
            target_radius_pixels = 1;
        end
        
        % Place target
        for dy = -target_radius_pixels:target_radius_pixels
            for dx = -target_radius_pixels:target_radius_pixels
                ix_target = ix_scene + dx;
                iy_target = iy_scene + dy;
                
                if ix_target > 0 && ix_target <= length(x_coords) && ...
                   iy_target > 0 && iy_target <= length(y_coords)
                    scene_matrix(iy_target, ix_target) = base_amplitude * target_intensity;
                end
            end
        end
    end
    
    % Create visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Main scene plot
    subplot(2, 3, 1);
    imagesc(x_coords, y_coords, scene_matrix);
    axis image;
    colormap gray;
    colorbar;
    xlabel('x (mm)');
    ylabel('y (mm)');
    title('Enhanced Target Scene (TheoryTest2-based)');
    set(gca, 'YDir', 'normal');
    
    % Target positions overlay
    subplot(2, 3, 2);
    imagesc(x_coords, y_coords, scene_matrix);
    axis image;
    colormap gray;
    hold on;
    for i = 1:size(target_positions, 1)
        plot(target_positions(i, 1), target_positions(i, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(target_positions(i, 1) + 5, target_positions(i, 2) + 5, sprintf('%d', i), 'Color', 'red', 'FontWeight', 'bold');
    end
    hold off;
    xlabel('x (mm)');
    ylabel('y (mm)');
    title('Target Positions (9 dots)');
    set(gca, 'YDir', 'normal');
    
    % Target size distribution
    subplot(2, 3, 3);
    target_sizes_used = target_sizes_mm(mod(0:8, length(target_sizes_mm)) + 1);
    bar(target_sizes_used);
    xlabel('Target Index');
    ylabel('Target Size (mm)');
    title('Target Size Distribution');
    grid on;
    
    % Target intensity distribution
    subplot(2, 3, 4);
    target_intensities_used = target_intensities(mod(0:8, length(target_intensities)) + 1);
    bar(target_intensities_used);
    xlabel('Target Index');
    ylabel('Target Intensity');
    title('Target Intensity Distribution');
    grid on;
    
    % Scene statistics
    subplot(2, 3, 5);
    scene_stats = struct();
    scene_stats.total_pixels = numel(scene_matrix);
    scene_stats.target_pixels = nnz(scene_matrix);
    scene_stats.background_pixels = scene_stats.total_pixels - scene_stats.target_pixels;
    scene_stats.target_coverage = scene_stats.target_pixels / scene_stats.total_pixels * 100;
    
    stats_data = [scene_stats.target_pixels, scene_stats.background_pixels];
    pie(stats_data, {'Target Pixels', 'Background Pixels'});
    title(sprintf('Scene Coverage (%.1f%% targets)', scene_stats.target_coverage));
    
    % Target spacing analysis
    subplot(2, 3, 6);
    target_distances = [];
    for i = 1:size(target_positions, 1)
        for j = i+1:size(target_positions, 1)
            dist = sqrt((target_positions(i,1) - target_positions(j,1))^2 + ...
                      (target_positions(i,2) - target_positions(j,2))^2);
            target_distances = [target_distances; dist];
        end
    end
    histogram(target_distances, 10);
    xlabel('Target Distance (mm)');
    ylabel('Frequency');
    title('Target Spacing Distribution');
    grid on;
    
    sgtitle('Enhanced Target Scene Analysis (TheoryTest2-based)', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'enhanced_target_scene.png'));
    close(gcf);
    
    % Save scene data
    save(fullfile(output_folder, 'enhanced_target_scene.mat'), 'scene_matrix', 'target_positions', 'scene_stats');
    
    fprintf('Enhanced target scene created with:\n');
    fprintf('  - 9 targets in 3x3 grid pattern\n');
    fprintf('  - Variable target sizes: %s mm\n', mat2str(target_sizes_mm));
    fprintf('  - Variable intensities: %s\n', mat2str(target_intensities));
    fprintf('  - Target coverage: %.1f%%\n', scene_stats.target_coverage);
    fprintf('  - Average target spacing: %.1f mm\n', mean(target_distances));
end

function visualize_sample_reconstructions(results, output_folder)
    % Visualize sample reconstructions from the sweep results
    
    if isempty(results) || ~isfield(results, 'configs')
        fprintf('No reconstruction results available for visualization\n');
        return;
    end
    
    % Find best and worst reconstructions
    psnr_values = zeros(length(results.metrics), 1);
    for i = 1:length(results.metrics)
        if ~isempty(results.metrics{i})
            psnr_values(i) = results.metrics{i}.psnr;
        end
    end
    
    [~, best_idx] = max(psnr_values);
    [~, worst_idx] = min(psnr_values);
    
    % Create visualization
    figure('Position', [100, 100, 1400, 1000]);
    
    % Best reconstruction
    subplot(2, 4, 1);
    if ~isempty(results.configs{best_idx})
        config = results.configs{best_idx};
        title(sprintf('Best Config (PSNR=%.2f dB)\nGrid: %.1fmm, Size: %.1fmm, Spacing: %.1fmm', ...
            psnr_values(best_idx), config.grid_step_mm, config.target_size_mm, config.grid_spacing_mm));
    else
        title('Best Configuration');
    end
    
    % Worst reconstruction
    subplot(2, 4, 2);
    if ~isempty(results.configs{worst_idx})
        config = results.configs{worst_idx};
        title(sprintf('Worst Config (PSNR=%.2f dB)\nGrid: %.1fmm, Size: %.1fmm, Spacing: %.1fmm', ...
            psnr_values(worst_idx), config.grid_step_mm, config.target_size_mm, config.grid_spacing_mm));
    else
        title('Worst Configuration');
    end
    
    % PSNR distribution
    subplot(2, 4, 3);
    histogram(psnr_values, 20);
    xlabel('PSNR (dB)');
    ylabel('Frequency');
    title('PSNR Distribution');
    grid on;
    
    % Parameter sensitivity plots
    subplot(2, 4, 4);
    grid_steps = zeros(length(results.configs), 1);
    target_sizes = zeros(length(results.configs), 1);
    for i = 1:length(results.configs)
        if ~isempty(results.configs{i})
            grid_steps(i) = results.configs{i}.grid_step_mm;
            target_sizes(i) = results.configs{i}.target_size_mm;
        end
    end
    scatter(grid_steps, target_sizes, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Step (mm)');
    ylabel('Target Size (mm)');
    title('PSNR vs Grid Parameters');
    grid on;
    
    % Target pattern analysis
    subplot(2, 4, 5);
    target_patterns = cell(length(results.configs), 1);
    for i = 1:length(results.configs)
        if ~isempty(results.configs{i})
            target_patterns{i} = results.configs{i}.target_pattern;
        end
    end
    % Handle empty or invalid target_patterns
    valid_patterns = target_patterns(~cellfun(@isempty, target_patterns));
    if ~isempty(valid_patterns)
        unique_patterns = unique(valid_patterns);
        pattern_psnr = zeros(length(unique_patterns), 1);
        for i = 1:length(unique_patterns)
            pattern_mask = strcmp(target_patterns, unique_patterns{i});
            pattern_psnr(i) = mean(psnr_values(pattern_mask));
        end
    else
        unique_patterns = {'unknown'};
        pattern_psnr = [mean(psnr_values)];
    end
    bar(pattern_psnr);
    set(gca, 'XTickLabel', unique_patterns);
    xlabel('Target Pattern');
    ylabel('Average PSNR (dB)');
    title('PSNR vs Target Pattern');
    grid on;
    
    % Correlation analysis
    subplot(2, 4, 6);
    correlation_values = zeros(length(results.metrics), 1);
    for i = 1:length(results.metrics)
        if ~isempty(results.metrics{i})
            correlation_values(i) = results.metrics{i}.correlation;
        end
    end
    scatter(psnr_values, correlation_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('PSNR (dB)');
    ylabel('Correlation');
    title('PSNR vs Correlation');
    grid on;
    
    % Time analysis
    subplot(2, 4, 7);
    reconstruction_times = zeros(length(results.times), 1);
    for i = 1:length(results.times)
        if ~isempty(results.times{i})
            reconstruction_times(i) = results.times{i}.total_time;
        end
    end
    scatter(reconstruction_times, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Reconstruction Time (s)');
    ylabel('PSNR (dB)');
    title('PSNR vs Time');
    grid on;
    
    % Performance summary
    subplot(2, 4, 8);
    performance_summary = struct();
    performance_summary.mean_psnr = mean(psnr_values);
    performance_summary.std_psnr = std(psnr_values);
    performance_summary.max_psnr = max(psnr_values);
    performance_summary.min_psnr = min(psnr_values);
    performance_summary.mean_correlation = mean(correlation_values);
    performance_summary.mean_time = mean(reconstruction_times);
    
    summary_text = sprintf(['Performance Summary:\n' ...
                           'Mean PSNR: %.2f ± %.2f dB\n' ...
                           'PSNR Range: %.2f - %.2f dB\n' ...
                           'Mean Correlation: %.4f\n' ...
                           'Mean Time: %.1f s'], ...
                           performance_summary.mean_psnr, performance_summary.std_psnr, ...
                           performance_summary.min_psnr, performance_summary.max_psnr, ...
                           performance_summary.mean_correlation, performance_summary.mean_time);
    
    text(0.1, 0.5, summary_text, 'FontSize', 10, 'VerticalAlignment', 'middle');
    axis off;
    
    sgtitle('Sample Reconstruction Analysis', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'sample_reconstructions.png'));
    close(gcf);
    
    % Save analysis data
    save(fullfile(output_folder, 'reconstruction_analysis.mat'), 'psnr_values', 'correlation_values', 'reconstruction_times', 'performance_summary');
end

function assess_target_quality(results, output_folder)
    % Assess target quality and reconstruction challenges
    
    fprintf('Assessing target quality and reconstruction challenges...\n');
    
    % Extract target-related parameters
    target_sizes = zeros(length(results.configs), 1);
    grid_spacings = zeros(length(results.configs), 1);
    grid_steps = zeros(length(results.configs), 1);
    psnr_values = zeros(length(results.configs), 1);
    
    for i = 1:length(results.configs)
        if ~isempty(results.configs{i})
            target_sizes(i) = results.configs{i}.target_size_mm;
            grid_spacings(i) = results.configs{i}.grid_spacing_mm;
            grid_steps(i) = results.configs{i}.grid_step_mm;
        end
        if ~isempty(results.metrics{i})
            psnr_values(i) = results.metrics{i}.psnr;
        end
    end
    
    % Calculate challenge metrics
    challenge_metrics = struct();
    
    % Target resolution challenge (smaller targets = more challenging)
    challenge_metrics.resolution_challenge = mean(target_sizes);
    
    % Target separation challenge (closer spacing = more challenging)
    challenge_metrics.separation_challenge = mean(grid_spacings);
    
    % Grid resolution challenge (finer grid = more challenging)
    challenge_metrics.grid_resolution_challenge = mean(grid_steps);
    
    % Overall challenge score (lower = more challenging)
    challenge_metrics.overall_challenge = (challenge_metrics.resolution_challenge + ...
                                         challenge_metrics.separation_challenge) / ...
                                         challenge_metrics.grid_resolution_challenge;
    
    % Create challenge assessment visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Challenge metrics
    subplot(2, 3, 1);
    challenge_data = [challenge_metrics.resolution_challenge, ...
                     challenge_metrics.separation_challenge, ...
                     challenge_metrics.grid_resolution_challenge];
    bar(challenge_data);
    set(gca, 'XTickLabel', {'Target Size', 'Target Spacing', 'Grid Step'});
    ylabel('Value (mm)');
    title('Challenge Metrics');
    grid on;
    
    % Target size vs PSNR
    subplot(2, 3, 2);
    scatter(target_sizes, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Target Size (mm)');
    ylabel('PSNR (dB)');
    title('Target Size vs PSNR');
    grid on;
    
    % Target spacing vs PSNR
    subplot(2, 3, 3);
    scatter(grid_spacings, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Target Spacing (mm)');
    ylabel('PSNR (dB)');
    title('Target Spacing vs PSNR');
    grid on;
    
    % Grid step vs PSNR
    subplot(2, 3, 4);
    scatter(grid_steps, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Step (mm)');
    ylabel('PSNR (dB)');
    title('Grid Step vs PSNR');
    grid on;
    
    % Challenge score distribution
    subplot(2, 3, 5);
    challenge_scores = (target_sizes + grid_spacings) ./ grid_steps;
    histogram(challenge_scores, 20);
    xlabel('Challenge Score');
    ylabel('Frequency');
    title('Challenge Score Distribution');
    grid on;
    
    % Challenge vs Performance
    subplot(2, 3, 6);
    scatter(challenge_scores, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Challenge Score');
    ylabel('PSNR (dB)');
    title('Challenge vs Performance');
    grid on;
    
    sgtitle('Target Quality Assessment', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'target_quality_assessment.png'));
    close(gcf);
    
    % Save challenge assessment
    save(fullfile(output_folder, 'challenge_assessment.mat'), 'challenge_metrics', 'challenge_scores');
    
    fprintf('Target quality assessment completed:\n');
    fprintf('  - Resolution challenge: %.2f mm\n', challenge_metrics.resolution_challenge);
    fprintf('  - Separation challenge: %.2f mm\n', challenge_metrics.separation_challenge);
    fprintf('  - Grid resolution challenge: %.2f mm\n', challenge_metrics.grid_resolution_challenge);
    fprintf('  - Overall challenge score: %.2f\n', challenge_metrics.overall_challenge);
end

function improved_targets = create_improved_target_configurations()
    % Create improved target configurations based on analysis
    
    improved_targets = struct();
    
    % Configuration 1: High Challenge (small targets, close spacing)
    improved_targets.high_challenge = struct();
    improved_targets.high_challenge.target_size_mm = 2;
    improved_targets.high_challenge.grid_spacing_mm = 10;
    improved_targets.high_challenge.grid_step_mm = 1;
    improved_targets.high_challenge.target_pattern = '3x3_grid';
    improved_targets.high_challenge.description = 'High challenge: 2mm targets, 10mm spacing, 1mm grid';
    
    % Configuration 2: Balanced Challenge
    improved_targets.balanced = struct();
    improved_targets.balanced.target_size_mm = 3;
    improved_targets.balanced.grid_spacing_mm = 15;
    improved_targets.balanced.grid_step_mm = 2;
    improved_targets.balanced.target_pattern = '3x3_grid';
    improved_targets.balanced.description = 'Balanced: 3mm targets, 15mm spacing, 2mm grid';
    
    % Configuration 3: Variable Challenge (different target sizes)
    improved_targets.variable = struct();
    improved_targets.variable.target_sizes_mm = [2, 3, 4, 5, 6, 7, 8, 9, 10];
    improved_targets.variable.grid_spacing_mm = 20;
    improved_targets.variable.grid_step_mm = 2;
    improved_targets.variable.target_pattern = '3x3_grid';
    improved_targets.variable.description = 'Variable: Different target sizes, 20mm spacing';
    
    % Configuration 4: Variable Intensity
    improved_targets.variable_intensity = struct();
    improved_targets.variable_intensity.target_size_mm = 4;
    improved_targets.variable_intensity.grid_spacing_mm = 18;
    improved_targets.variable_intensity.grid_step_mm = 2;
    improved_targets.variable_intensity.target_intensities = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3];
    improved_targets.variable_intensity.target_pattern = '3x3_grid';
    improved_targets.variable_intensity.description = 'Variable intensity: 4mm targets with different amplitudes';
    
    % Configuration 5: Optimal from Analysis
    improved_targets.optimal = struct();
    improved_targets.optimal.target_size_mm = 4;
    improved_targets.optimal.grid_spacing_mm = 20;
    improved_targets.optimal.grid_step_mm = 2;
    improved_targets.optimal.target_pattern = '3x3_grid';
    improved_targets.optimal.description = 'Optimal: Based on sweep analysis results';
    
    fprintf('Created %d improved target configurations\n', length(fieldnames(improved_targets)));
end

function create_comprehensive_dashboard(results, output_folder)
    % Create comprehensive visualization dashboard
    
    if isempty(results) || ~isfield(results, 'configs')
        fprintf('No results available for dashboard\n');
        return;
    end
    
    % Extract all metrics
    psnr_values = zeros(length(results.metrics), 1);
    correlation_values = zeros(length(results.metrics), 1);
    reconstruction_times = zeros(length(results.times), 1);
    
    grid_steps = zeros(length(results.configs), 1);
    target_sizes = zeros(length(results.configs), 1);
    grid_spacings = zeros(length(results.configs), 1);
    num_acquisitions = zeros(length(results.configs), 1);
    lambda_tv_values = zeros(length(results.configs), 1);
    
    for i = 1:length(results.metrics)
        if ~isempty(results.metrics{i})
            psnr_values(i) = results.metrics{i}.psnr;
            correlation_values(i) = results.metrics{i}.correlation;
        end
        if ~isempty(results.times{i})
            reconstruction_times(i) = results.times{i}.total_time;
        end
        if ~isempty(results.configs{i})
            config = results.configs{i};
            grid_steps(i) = config.grid_step_mm;
            target_sizes(i) = config.target_size_mm;
            grid_spacings(i) = config.grid_spacing_mm;
            num_acquisitions(i) = config.num_acquisitions;
            lambda_tv_values(i) = config.lambda_tv_reg;
        end
    end
    
    % Create comprehensive dashboard
    figure('Position', [50, 50, 1600, 1200]);
    
    % Performance overview
    subplot(3, 4, 1);
    performance_data = [mean(psnr_values), mean(correlation_values), mean(reconstruction_times)];
    bar(performance_data);
    set(gca, 'XTickLabel', {'PSNR (dB)', 'Correlation', 'Time (s)'});
    title('Average Performance Metrics');
    grid on;
    
    % PSNR distribution
    subplot(3, 4, 2);
    histogram(psnr_values, 20, 'Normalization', 'probability');
    xlabel('PSNR (dB)');
    ylabel('Probability');
    title('PSNR Distribution');
    grid on;
    
    % Correlation distribution
    subplot(3, 4, 3);
    histogram(correlation_values, 20, 'Normalization', 'probability');
    xlabel('Correlation');
    ylabel('Probability');
    title('Correlation Distribution');
    grid on;
    
    % Time distribution
    subplot(3, 4, 4);
    histogram(reconstruction_times, 20, 'Normalization', 'probability');
    xlabel('Reconstruction Time (s)');
    ylabel('Probability');
    title('Time Distribution');
    grid on;
    
    % Parameter sensitivity plots
    subplot(3, 4, 5);
    scatter(grid_steps, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Step (mm)');
    ylabel('PSNR (dB)');
    title('Grid Step vs PSNR');
    grid on;
    
    subplot(3, 4, 6);
    scatter(target_sizes, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Target Size (mm)');
    ylabel('PSNR (dB)');
    title('Target Size vs PSNR');
    grid on;
    
    subplot(3, 4, 7);
    scatter(grid_spacings, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Grid Spacing (mm)');
    ylabel('PSNR (dB)');
    title('Grid Spacing vs PSNR');
    grid on;
    
    subplot(3, 4, 8);
    scatter(num_acquisitions, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Number of Acquisitions');
    ylabel('PSNR (dB)');
    title('Acquisitions vs PSNR');
    grid on;
    
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
    
    % Best configuration summary
    subplot(3, 4, 11);
    [~, best_idx] = max(psnr_values);
    if ~isempty(results.configs{best_idx})
        best_config = results.configs{best_idx};
        best_params = [best_config.grid_step_mm, best_config.target_size_mm, ...
                      best_config.grid_spacing_mm, best_config.num_acquisitions, ...
                      best_config.lambda_tv_reg];
        param_names = {'Grid Step', 'Target Size', 'Grid Spacing', 'Acquisitions', 'Lambda TV'};
        bar(best_params);
        set(gca, 'XTickLabel', param_names);
        ylabel('Value');
        title(sprintf('Best Configuration\n(PSNR=%.2f dB)', psnr_values(best_idx)));
        grid on;
    end
    
    % Challenge assessment
    subplot(3, 4, 12);
    challenge_scores = (target_sizes + grid_spacings) ./ grid_steps;
    scatter(challenge_scores, psnr_values, 50, psnr_values, 'filled');
    colorbar;
    xlabel('Challenge Score');
    ylabel('PSNR (dB)');
    title('Challenge vs Performance');
    grid on;
    
    sgtitle('Comprehensive Analysis Dashboard', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'comprehensive_dashboard.png'));
    close(gcf);
    
    fprintf('Comprehensive dashboard created\n');
end

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
            
            param_combinations = [param_combinations; config];
        end
    end
    
    fprintf('Reconstructed results for %d runs\n', num_runs);
end 