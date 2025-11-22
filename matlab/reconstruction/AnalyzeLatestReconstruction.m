% Analyze Latest Reconstruction Results
% This script loads and analyzes the latest reconstruction to see if it's actually good
clearvars;
clc;
close all;

%% Load the latest reconstruction results
fprintf('--- Analyzing Latest Reconstruction Results ---\n');

% Load reconstruction results (updated to latest folder)
reconstruction_file = '/Users/deepshikhakaul/Documents/MATLAB/Optimized_071925_012/optimized_reconstruction_results.mat';
if ~exist(reconstruction_file, 'file')
    error('Reconstruction file not found: %s', reconstruction_file);
end

fprintf('Loading reconstruction results from: %s\n', reconstruction_file);
reconstruction_data = load(reconstruction_file);

% Extract reconstruction results
if isfield(reconstruction_data, 'image_recon')
    image_recon = reconstruction_data.image_recon;
    fprintf('Reconstruction image loaded: %s\n', mat2str(size(image_recon)));
else
    error('No image_recon field found in reconstruction file');
end

if isfield(reconstruction_data, 'convergence_history')
    convergence_history = reconstruction_data.convergence_history;
    fprintf('Convergence history loaded: %d iterations\n', length(convergence_history));
else
    convergence_history = [];
    fprintf('No convergence history found\n');
end

if isfield(reconstruction_data, 'params')
    params = reconstruction_data.params;
    fprintf('Parameters loaded\n');
else
    params = [];
    fprintf('No parameters found\n');
end

%% Analyze reconstruction quality
fprintf('\n--- Reconstruction Quality Analysis ---\n');

% Basic statistics
image_stats = struct();
image_stats.mean = mean(image_recon(:));
image_stats.std = std(image_recon(:));
image_stats.min = min(image_recon(:));
image_stats.max = max(image_recon(:));
image_stats.range = image_stats.max - image_stats.min;
image_stats.dynamic_range_db = 20 * log10(image_stats.max / (abs(image_stats.min) + eps));

fprintf('Image Statistics:\n');
fprintf('  Mean: %.6f\n', image_stats.mean);
fprintf('  Std: %.6f\n', image_stats.std);
fprintf('  Min: %.6f\n', image_stats.min);
fprintf('  Max: %.6f\n', image_stats.max);
fprintf('  Range: %.6f\n', image_stats.range);
fprintf('  Dynamic Range: %.2f dB\n', image_stats.dynamic_range_db);

% Check for meaningful features
nonzero_pixels = sum(image_recon(:) ~= 0);
total_pixels = numel(image_recon);
sparsity = 1 - (nonzero_pixels / total_pixels);

fprintf('\nFeature Analysis:\n');
fprintf('  Total pixels: %d\n', total_pixels);
fprintf('  Non-zero pixels: %d\n', nonzero_pixels);
fprintf('  Sparsity: %.2f%%\n', sparsity * 100);

% Check for gradient-like behavior
[gradient_x, gradient_y] = gradient(image_recon);
gradient_magnitude = sqrt(gradient_x.^2 + gradient_y.^2);
avg_gradient = mean(gradient_magnitude(:));

fprintf('  Average gradient magnitude: %.6f\n', avg_gradient);

% Check for distinct features
[peaks, peak_locs] = findpeaks(image_recon(:), 'MinPeakHeight', image_stats.max * 0.1);
num_peaks = length(peaks);
fprintf('  Number of significant peaks: %d\n', num_peaks);

%% Visualize the reconstruction
figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Main reconstruction image
subplot(2, 3, 1);
imagesc(image_recon);
colorbar;
title('Reconstructed Image');
xlabel('X (pixels)'); ylabel('Z (pixels)');
axis equal;

% Histogram of pixel values
subplot(2, 3, 2);
histogram(image_recon(:), 50);
title('Pixel Value Distribution');
xlabel('Pixel Value'); ylabel('Count');
grid on;

% Line plot through center
subplot(2, 3, 3);
center_row = round(size(image_recon, 1) / 2);
plot(image_recon(center_row, :), 'b-', 'LineWidth', 2);
title(sprintf('Center Row (Row %d)', center_row));
xlabel('X (pixels)'); ylabel('Intensity');
grid on;

% Convergence history
if ~isempty(convergence_history)
    subplot(2, 3, 4);
    semilogy(convergence_history, 'r-', 'LineWidth', 2);
    title('ADMM Convergence History');
    xlabel('Iteration'); ylabel('Residual');
    grid on;
end

% 3D surface plot
subplot(2, 3, 5);
surf(image_recon);
title('3D Surface Plot');
xlabel('X'); ylabel('Z'); zlabel('Intensity');
shading interp;

% Gradient magnitude
subplot(2, 3, 6);
imagesc(gradient_magnitude);
colorbar;
title('Gradient Magnitude');
xlabel('X (pixels)'); ylabel('Z (pixels)');
axis equal;

sgtitle('Latest Reconstruction Analysis (Fixed Time Delays)');
set(gcf, 'Color', 'w');

%% Compare with previous results
fprintf('\n--- Comparison with Previous Results ---\n');

% Try to load a previous result for comparison
previous_file = '/Users/deepshikhakaul/Documents/MATLAB/Optimized_071925_011/optimized_reconstruction_results.mat';
if exist(previous_file, 'file')
    fprintf('Loading previous results for comparison...\n');
    previous_data = load(previous_file);
    
    if isfield(previous_data, 'image_recon')
        previous_image = previous_data.image_recon;
        
        % Calculate improvement metrics
        current_energy = sum(image_recon(:).^2);
        previous_energy = sum(previous_image(:).^2);
        energy_ratio = current_energy / previous_energy;
        
        current_contrast = (max(image_recon(:)) - min(image_recon(:))) / (max(image_recon(:)) + min(image_recon(:)) + eps);
        previous_contrast = (max(previous_image(:)) - min(previous_image(:))) / (max(previous_image(:)) + min(previous_image(:)) + eps);
        contrast_ratio = current_contrast / previous_contrast;
        
        fprintf('Improvement Analysis:\n');
        fprintf('  Energy ratio (current/previous): %.2f\n', energy_ratio);
        fprintf('  Contrast ratio (current/previous): %.2f\n', contrast_ratio);
        
        if energy_ratio > 1
            fprintf('  ✓ Current reconstruction has more energy\n');
        else
            fprintf('  ✗ Previous reconstruction had more energy\n');
        end
        
        if contrast_ratio > 1
            fprintf('  ✓ Current reconstruction has better contrast\n');
        else
            fprintf('  ✗ Previous reconstruction had better contrast\n');
        end
    end
else
    fprintf('No previous results found for comparison\n');
end

%% Summary assessment
fprintf('\n--- Summary Assessment ---\n');

% Assess reconstruction quality
quality_score = 0;
max_score = 5;

% Check dynamic range
if image_stats.dynamic_range_db > 20
    quality_score = quality_score + 1;
    fprintf('✓ Good dynamic range (%.1f dB)\n', image_stats.dynamic_range_db);
else
    fprintf('✗ Poor dynamic range (%.1f dB)\n', image_stats.dynamic_range_db);
end

% Check for features (not just gradient)
if num_peaks > 5
    quality_score = quality_score + 1;
    fprintf('✓ Multiple distinct features (%d peaks)\n', num_peaks);
else
    fprintf('✗ Few distinct features (%d peaks)\n', num_peaks);
end

% Check sparsity (should have some structure)
if sparsity < 0.9
    quality_score = quality_score + 1;
    fprintf('✓ Good feature density (%.1f%% non-zero)\n', (1-sparsity)*100);
else
    fprintf('✗ Too sparse (%.1f%% non-zero)\n', (1-sparsity)*100);
end

% Check convergence
if ~isempty(convergence_history) && length(convergence_history) < 50
    quality_score = quality_score + 1;
    fprintf('✓ Good convergence (%d iterations)\n', length(convergence_history));
else
    fprintf('✗ Poor convergence or too many iterations\n');
end

% Check for reasonable values
if image_stats.std > 0.001
    quality_score = quality_score + 1;
    fprintf('✓ Reasonable pixel variation (std=%.6f)\n', image_stats.std);
else
    fprintf('✗ Low pixel variation (std=%.6f)\n', image_stats.std);
end

fprintf('\nOverall Quality Score: %d/%d\n', quality_score, max_score);

if quality_score >= 4
    fprintf('🎉 EXCELLENT reconstruction quality!\n');
elseif quality_score >= 3
    fprintf('✅ GOOD reconstruction quality\n');
elseif quality_score >= 2
    fprintf('⚠️  FAIR reconstruction quality\n');
else
    fprintf('❌ POOR reconstruction quality\n');
end

fprintf('\n--- Analysis Complete ---\n'); 