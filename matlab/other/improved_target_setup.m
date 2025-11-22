%% ===== IMPROVED TARGET SETUP BASED ON THEORYTEST2 =====
% Enhanced target configuration with better challenge and visualization
% Designed for optimal reconstruction testing with 9 small, clearly distinguishable targets

clear; clc; close all;

%% ===== CONFIGURATION =====
% Output setup
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('improved_target_setup', timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fprintf('=== IMPROVED TARGET SETUP ===\n');
fprintf('Based on TheoryTest2 approach with enhanced challenge\n');
fprintf('Saving results to: %s\n\n', output_folder);

%% ===== ENHANCED TARGET CONFIGURATIONS =====
fprintf('Creating enhanced target configurations...\n');

% Configuration 1: High Challenge (TheoryTest2-inspired)
config_high_challenge = create_high_challenge_config();
visualize_target_config(config_high_challenge, 'High Challenge', output_folder);

% Configuration 2: Variable Challenge
config_variable = create_variable_challenge_config();
visualize_target_config(config_variable, 'Variable Challenge', output_folder);

% Configuration 3: Optimal Challenge
config_optimal = create_optimal_challenge_config();
visualize_target_config(config_optimal, 'Optimal Challenge', output_folder);

% Configuration 4: Realistic Challenge
config_realistic = create_realistic_challenge_config();
visualize_target_config(config_realistic, 'Realistic Challenge', output_folder);

%% ===== COMPARATIVE ANALYSIS =====
fprintf('Creating comparative analysis...\n');

configs = {config_high_challenge, config_variable, config_optimal, config_realistic};
config_names = {'High Challenge', 'Variable Challenge', 'Optimal Challenge', 'Realistic Challenge'};

create_comparative_analysis(configs, config_names, output_folder);

%% ===== RECOMMENDATIONS =====
fprintf('Generating recommendations...\n');

create_recommendations(configs, config_names, output_folder);

fprintf('\n=== IMPROVED TARGET SETUP COMPLETE ===\n');
fprintf('Results saved to: %s\n', output_folder);

%% ===== HELPER FUNCTIONS =====

function config = create_high_challenge_config()
    % Create high challenge configuration based on TheoryTest2
    
    config = struct();
    
    % Grid parameters (TheoryTest2-inspired)
    config.grid_width_mm = 150;
    config.grid_step_mm = 1;  % Very fine resolution
    config.target_distance_mm = 250;
    
    % Target positions (3x3 grid like TheoryTest2)
    config.target_positions = [
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
    
    % High challenge parameters
    config.target_sizes_mm = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];  % Very small targets
    config.target_intensities = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4];  % Variable intensity
    config.target_spacing_mm = 15;  % Close spacing
    
    config.description = 'High Challenge: Very small targets (1-5mm), close spacing (15mm), fine grid (1mm)';
    config.challenge_level = 'High';
    config.expected_difficulty = 'Very challenging reconstruction';
end

function config = create_variable_challenge_config()
    % Create variable challenge configuration
    
    config = struct();
    
    % Grid parameters
    config.grid_width_mm = 150;
    config.grid_step_mm = 2;
    config.target_distance_mm = 250;
    
    % Target positions (same 3x3 grid)
    config.target_positions = [
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
    
    % Variable challenge parameters
    config.target_sizes_mm = [2, 3, 4, 5, 6, 7, 8, 9, 10];  % Wide range of sizes
    config.target_intensities = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3];  % Variable intensity
    config.target_spacing_mm = 20;  % Moderate spacing
    
    config.description = 'Variable Challenge: Different target sizes (2-10mm), variable intensity, moderate spacing (20mm)';
    config.challenge_level = 'Variable';
    config.expected_difficulty = 'Mixed difficulty based on target properties';
end

function config = create_optimal_challenge_config()
    % Create optimal challenge configuration based on analysis
    
    config = struct();
    
    % Grid parameters
    config.grid_width_mm = 150;
    config.grid_step_mm = 2;
    config.target_distance_mm = 250;
    
    % Target positions (same 3x3 grid)
    config.target_positions = [
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
    
    % Optimal challenge parameters (based on sweep analysis)
    config.target_sizes_mm = [3, 3, 4, 4, 4, 5, 5, 5, 6];  % Balanced sizes
    config.target_intensities = [0.8, 0.9, 1.0, 1.0, 1.1, 1.1, 1.2, 1.2, 1.3];  % Balanced intensity
    config.target_spacing_mm = 18;  % Optimal spacing
    
    config.description = 'Optimal Challenge: Balanced target sizes (3-6mm), balanced intensity, optimal spacing (18mm)';
    config.challenge_level = 'Optimal';
    config.expected_difficulty = 'Balanced challenge for good reconstruction';
end

function config = create_realistic_challenge_config()
    % Create realistic challenge configuration
    
    config = struct();
    
    % Grid parameters
    config.grid_width_mm = 150;
    config.grid_step_mm = 2;
    config.target_distance_mm = 250;
    
    % Target positions (same 3x3 grid)
    config.target_positions = [
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
    
    % Realistic challenge parameters
    config.target_sizes_mm = [4, 4, 4, 4, 4, 4, 4, 4, 4];  % Uniform size
    config.target_intensities = [0.9, 0.9, 0.9, 0.9, 1.0, 0.9, 0.9, 0.9, 0.9];  % Nearly uniform
    config.target_spacing_mm = 20;  % Standard spacing
    
    config.description = 'Realistic Challenge: Uniform target sizes (4mm), nearly uniform intensity, standard spacing (20mm)';
    config.challenge_level = 'Realistic';
    config.expected_difficulty = 'Realistic challenge for practical applications';
end

function visualize_target_config(config, config_name, output_folder)
    % Visualize target configuration
    
    % Create imaging grid
    x_coords = -config.grid_width_mm/2 : config.grid_step_mm : config.grid_width_mm/2;
    y_coords = -config.grid_width_mm/2 : config.grid_step_mm : config.grid_width_mm/2;
    [X_mesh, Y_mesh] = meshgrid(x_coords, y_coords);
    
    % Create scene matrix
    scene_matrix = zeros(length(y_coords), length(x_coords));
    
    % Place targets
    for i = 1:size(config.target_positions, 1)
        x_pos_mm = config.target_positions(i, 1);
        y_pos_mm = config.target_positions(i, 2);
        base_amplitude = config.target_positions(i, 3);
        
        % Get target properties
        target_size_mm = config.target_sizes_mm(i);
        target_intensity = config.target_intensities(i);
        
        % Find nearest grid point
        [~, ix_scene] = min(abs(x_coords - x_pos_mm));
        [~, iy_scene] = min(abs(y_coords - y_pos_mm));
        
        % Create target with specified size
        target_radius_pixels = round(target_size_mm / (2 * config.grid_step_mm));
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
    title(sprintf('%s Target Scene', config_name));
    set(gca, 'YDir', 'normal');
    
    % Target positions overlay
    subplot(2, 3, 2);
    imagesc(x_coords, y_coords, scene_matrix);
    axis image;
    colormap gray;
    hold on;
    for i = 1:size(config.target_positions, 1)
        plot(config.target_positions(i, 1), config.target_positions(i, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(config.target_positions(i, 1) + 3, config.target_positions(i, 2) + 3, sprintf('%d', i), 'Color', 'red', 'FontWeight', 'bold');
    end
    hold off;
    xlabel('x (mm)');
    ylabel('y (mm)');
    title('Target Positions');
    set(gca, 'YDir', 'normal');
    
    % Target size distribution
    subplot(2, 3, 3);
    bar(config.target_sizes_mm);
    xlabel('Target Index');
    ylabel('Target Size (mm)');
    title('Target Size Distribution');
    grid on;
    
    % Target intensity distribution
    subplot(2, 3, 4);
    bar(config.target_intensities);
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
    
    % Configuration summary
    subplot(2, 3, 6);
    summary_text = sprintf(['Configuration Summary:\n' ...
                           'Challenge Level: %s\n' ...
                           'Grid Step: %.1f mm\n' ...
                           'Target Spacing: %.1f mm\n' ...
                           'Avg Target Size: %.1f mm\n' ...
                           'Avg Intensity: %.2f\n' ...
                           'Expected Difficulty: %s'], ...
                           config.challenge_level, config.grid_step_mm, ...
                           config.target_spacing_mm, mean(config.target_sizes_mm), ...
                           mean(config.target_intensities), config.expected_difficulty);
    
    text(0.1, 0.5, summary_text, 'FontSize', 10, 'VerticalAlignment', 'middle');
    axis off;
    
    sgtitle(sprintf('%s Configuration', config_name), 'FontSize', 16);
    
    % Save figure
    safe_filename = strrep(config_name, ' ', '_');
    saveas(gcf, fullfile(output_folder, sprintf('%s_target_config.png', safe_filename)));
    close(gcf);
    
    % Save configuration data
    save(fullfile(output_folder, sprintf('%s_config.mat', safe_filename)), 'config', 'scene_matrix', 'scene_stats');
    
    fprintf('Created %s configuration:\n', config_name);
    fprintf('  - %s\n', config.description);
    fprintf('  - Target coverage: %.1f%%\n', scene_stats.target_coverage);
    fprintf('  - Average target size: %.1f mm\n', mean(config.target_sizes_mm));
    fprintf('  - Average intensity: %.2f\n', mean(config.target_intensities));
end

function create_comparative_analysis(configs, config_names, output_folder)
    % Create comparative analysis of all configurations
    
    % Extract metrics for comparison
    metrics = struct();
    
    for i = 1:length(configs)
        config = configs{i};
        
        % Calculate metrics
        metrics(i).name = config_names{i};
        metrics(i).challenge_level = config.challenge_level;
        metrics(i).avg_target_size = mean(config.target_sizes_mm);
        metrics(i).std_target_size = std(config.target_sizes_mm);
        metrics(i).avg_intensity = mean(config.target_intensities);
        metrics(i).std_intensity = std(config.target_intensities);
        metrics(i).target_spacing = config.target_spacing_mm;
        metrics(i).grid_step = config.grid_step_mm;
        
        % Calculate challenge score (lower = more challenging)
        metrics(i).challenge_score = (metrics(i).avg_target_size + metrics(i).target_spacing) / metrics(i).grid_step;
    end
    
    % Create comparative visualization
    figure('Position', [100, 100, 1400, 1000]);
    
    % Challenge score comparison
    subplot(2, 3, 1);
    challenge_scores = [metrics.challenge_score];
    bar(challenge_scores);
    set(gca, 'XTickLabel', config_names);
    ylabel('Challenge Score');
    title('Challenge Score Comparison');
    grid on;
    
    % Target size comparison
    subplot(2, 3, 2);
    avg_sizes = [metrics.avg_target_size];
    std_sizes = [metrics.std_target_size];
    bar(avg_sizes);
    hold on;
    errorbar(1:length(avg_sizes), avg_sizes, std_sizes, 'k.', 'LineWidth', 2);
    hold off;
    set(gca, 'XTickLabel', config_names);
    ylabel('Average Target Size (mm)');
    title('Target Size Comparison');
    grid on;
    
    % Intensity comparison
    subplot(2, 3, 3);
    avg_intensities = [metrics.avg_intensity];
    std_intensities = [metrics.std_intensity];
    bar(avg_intensities);
    hold on;
    errorbar(1:length(avg_intensities), avg_intensities, std_intensities, 'k.', 'LineWidth', 2);
    hold off;
    set(gca, 'XTickLabel', config_names);
    ylabel('Average Intensity');
    title('Intensity Comparison');
    grid on;
    
    % Spacing comparison
    subplot(2, 3, 4);
    spacings = [metrics.target_spacing];
    bar(spacings);
    set(gca, 'XTickLabel', config_names);
    ylabel('Target Spacing (mm)');
    title('Target Spacing Comparison');
    grid on;
    
    % Grid resolution comparison
    subplot(2, 3, 5);
    grid_steps = [metrics.grid_step];
    bar(grid_steps);
    set(gca, 'XTickLabel', config_names);
    ylabel('Grid Step (mm)');
    title('Grid Resolution Comparison');
    grid on;
    
    % Challenge level summary
    subplot(2, 3, 6);
    challenge_levels = {metrics.challenge_level};
    level_counts = zeros(1, 4);
    unique_levels = unique(challenge_levels);
    for i = 1:length(unique_levels)
        level_counts(i) = sum(strcmp(challenge_levels, unique_levels{i}));
    end
    pie(level_counts, unique_levels);
    title('Challenge Level Distribution');
    
    sgtitle('Comparative Analysis of Target Configurations', 'FontSize', 16);
    
    saveas(gcf, fullfile(output_folder, 'comparative_analysis.png'));
    close(gcf);
    
    % Save comparative data
    save(fullfile(output_folder, 'comparative_analysis.mat'), 'metrics');
    
    fprintf('Comparative analysis completed:\n');
    for i = 1:length(metrics)
        fprintf('  %s: Challenge Score=%.2f, Avg Size=%.1f±%.1f mm, Avg Intensity=%.2f±%.2f\n', ...
            metrics(i).name, metrics(i).challenge_score, metrics(i).avg_target_size, metrics(i).std_target_size, ...
            metrics(i).avg_intensity, metrics(i).std_intensity);
    end
end

function create_recommendations(configs, config_names, output_folder)
    % Create recommendations based on analysis
    
    fprintf('Generating recommendations...\n');
    
    % Create recommendations file
    recommendations_file = fullfile(output_folder, 'recommendations.txt');
    fid = fopen(recommendations_file, 'w');
    
    fprintf(fid, '=== IMPROVED TARGET SETUP RECOMMENDATIONS ===\n\n');
    fprintf(fid, 'Based on TheoryTest2 approach and parameter sweep analysis\n\n');
    
    % Configuration recommendations
    fprintf(fid, 'CONFIGURATION RECOMMENDATIONS:\n\n');
    
    for i = 1:length(configs)
        config = configs{i};
        fprintf(fid, '%d. %s Configuration:\n', i, config_names{i});
        fprintf(fid, '   - Challenge Level: %s\n', config.challenge_level);
        fprintf(fid, '   - Description: %s\n', config.description);
        fprintf(fid, '   - Target Sizes: %s mm\n', mat2str(config.target_sizes_mm));
        fprintf(fid, '   - Target Intensities: %s\n', mat2str(config.target_intensities));
        fprintf(fid, '   - Grid Step: %.1f mm\n', config.grid_step_mm);
        fprintf(fid, '   - Target Spacing: %.1f mm\n', config.target_spacing_mm);
        fprintf(fid, '   - Expected Difficulty: %s\n\n', config.expected_difficulty);
    end
    
    % Usage recommendations
    fprintf(fid, 'USAGE RECOMMENDATIONS:\n\n');
    fprintf(fid, '1. High Challenge Configuration:\n');
    fprintf(fid, '   - Use for: Testing reconstruction algorithm limits\n');
    fprintf(fid, '   - Best for: Algorithm development and optimization\n');
    fprintf(fid, '   - Expected: Very challenging reconstruction\n\n');
    
    fprintf(fid, '2. Variable Challenge Configuration:\n');
    fprintf(fid, '   - Use for: Testing algorithm robustness\n');
    fprintf(fid, '   - Best for: Comprehensive algorithm evaluation\n');
    fprintf(fid, '   - Expected: Mixed difficulty results\n\n');
    
    fprintf(fid, '3. Optimal Challenge Configuration:\n');
    fprintf(fid, '   - Use for: Standard performance testing\n');
    fprintf(fid, '   - Best for: Balanced challenge and performance\n');
    fprintf(fid, '   - Expected: Good reconstruction quality\n\n');
    
    fprintf(fid, '4. Realistic Challenge Configuration:\n');
    fprintf(fid, '   - Use for: Practical application testing\n');
    fprintf(fid, '   - Best for: Real-world scenario simulation\n');
    fprintf(fid, '   - Expected: Realistic performance assessment\n\n');
    
    % Implementation recommendations
    fprintf(fid, 'IMPLEMENTATION RECOMMENDATIONS:\n\n');
    fprintf(fid, '1. Grid Setup:\n');
    fprintf(fid, '   - Use 150mm x 150mm imaging grid\n');
    fprintf(fid, '   - Position targets at 250mm distance\n');
    fprintf(fid, '   - Use 3x3 grid pattern for 9 targets\n\n');
    
    fprintf(fid, '2. Target Properties:\n');
    fprintf(fid, '   - Vary target sizes for challenge (1-10mm range)\n');
    fprintf(fid, '   - Use variable intensities (0.5-1.4 range)\n');
    fprintf(fid, '   - Maintain clear background (zero amplitude)\n\n');
    
    fprintf(fid, '3. Spacing Considerations:\n');
    fprintf(fid, '   - Close spacing (15-20mm) for high challenge\n');
    fprintf(fid, '   - Moderate spacing (20-25mm) for balanced challenge\n');
    fprintf(fid, '   - Ensure targets are clearly distinguishable\n\n');
    
    fprintf(fid, '4. Grid Resolution:\n');
    fprintf(fid, '   - Fine grid (1-2mm) for high resolution\n');
    fprintf(fid, '   - Coarse grid (3-5mm) for computational efficiency\n');
    fprintf(fid, '   - Balance resolution with computation time\n\n');
    
    % Best practices
    fprintf(fid, 'BEST PRACTICES:\n\n');
    fprintf(fid, '1. Always visualize targets before reconstruction\n');
    fprintf(fid, '2. Test multiple configurations for robustness\n');
    fprintf(fid, '3. Monitor target coverage percentage\n');
    fprintf(fid, '4. Ensure targets are within imaging grid bounds\n');
    fprintf(fid, '5. Use consistent target patterns for comparison\n');
    fprintf(fid, '6. Document configuration parameters for reproducibility\n\n');
    
    fclose(fid);
    
    fprintf('Recommendations saved to: %s\n', recommendations_file);
end 