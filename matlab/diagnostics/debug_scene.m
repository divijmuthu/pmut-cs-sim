%% debug_scene.m - Debug the scene creation
clear; clc; close all;

% Configuration
config = struct();
config.grid_step_m = 0.0075;       % 7.5 mm grid step
config.target_size_mm = 3.0;       % 3.0 mm targets
config.grid_spacing_mm = 20.0;     % 20.0 mm spacing
config.grid_width_m = 0.150;       % Imaging grid width (m)
config.target_distance_m = 0.150;  % Target distance (m)
config.grid_depth_range_m = 0.100; % Grid depth range (m)

% Create imaging grid
x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);

fprintf('Grid size: %d x %d = %d pixels\n', size(X_mesh, 1), size(X_mesh, 2), numel(X_mesh));
fprintf('X range: %.1f to %.1f mm\n', min(x_coords_img)*1000, max(x_coords_img)*1000);
fprintf('Z range: %.1f to %.1f mm\n', min(z_coords_img)*1000, max(z_coords_img)*1000);

% Create scene
scene_matrix = zeros(size(X_mesh));

% Convert to mm
X_mm = X_mesh * 1000;
Z_mm = Z_mesh * 1000;

% Grid parameters
target_size_mm = config.target_size_mm;
grid_spacing_mm = config.grid_spacing_mm;
grid_step_mm = config.grid_step_m * 1000;
target_radius_pixels = round(target_size_mm / (2 * grid_step_mm));

fprintf('Grid step: %.1f mm\n', grid_step_mm);
fprintf('Target radius in pixels: %d\n', target_radius_pixels);

if target_radius_pixels < 1
    target_radius_pixels = 1;
    fprintf('Adjusted target radius to: %d pixels\n', target_radius_pixels);
end

% 3x3 grid positions
grid_offset_x_mm = 0;
grid_offset_z_mm = 150;

target_positions = [];
for row = 1:3
    for col = 1:3
        x_pos_mm = grid_offset_x_mm + (col - 2) * grid_spacing_mm;
        z_pos_mm = grid_offset_z_mm + (row - 2) * grid_spacing_mm;
        target_positions = [target_positions; x_pos_mm, z_pos_mm, 1.0];
    end
end

fprintf('\nTarget positions:\n');
for i = 1:size(target_positions, 1)
    fprintf('Target %d: (%.1f, %.1f) mm\n', i, target_positions(i, 1), target_positions(i, 2));
end

% Place targets
for i = 1:size(target_positions, 1)
    x_pos_mm = target_positions(i, 1);
    z_pos_mm = target_positions(i, 2);
    amplitude = target_positions(i, 3);
    
    [~, ix_scene] = min(abs(X_mm(1,:) - x_pos_mm));
    [~, iz_scene] = min(abs(Z_mm(:,1) - z_pos_mm));
    
    fprintf('Target %d: Placing at grid indices (%d, %d)\n', i, iz_scene, ix_scene);
    
    for dz = -target_radius_pixels:target_radius_pixels
        for dx = -target_radius_pixels:target_radius_pixels
            if iz_scene + dz >= 1 && iz_scene + dz <= size(scene_matrix, 1) && ...
               ix_scene + dx >= 1 && ix_scene + dx <= size(scene_matrix, 2)
                scene_matrix(iz_scene + dz, ix_scene + dx) = amplitude;
            end
        end
    end
end

% Analyze scene matrix
fprintf('\nScene matrix analysis:\n');
fprintf('Scene matrix size: %d x %d\n', size(scene_matrix, 1), size(scene_matrix, 2));
fprintf('Min value: %.6f\n', min(scene_matrix(:)));
fprintf('Max value: %.6f\n', max(scene_matrix(:)));
fprintf('Number of non-zero pixels: %d\n', nnz(scene_matrix));
fprintf('Total pixels: %d\n', numel(scene_matrix));
fprintf('Non-zero percentage: %.2f%%\n', nnz(scene_matrix)/numel(scene_matrix)*100);

% Show where the non-zero pixels are
[nonzero_rows, nonzero_cols] = find(scene_matrix > 0);
if ~isempty(nonzero_rows)
    fprintf('\nNon-zero pixel locations:\n');
    for i = 1:min(10, length(nonzero_rows))
        fprintf('Pixel (%d, %d) = %.1f\n', nonzero_rows(i), nonzero_cols(i), scene_matrix(nonzero_rows(i), nonzero_cols(i)));
    end
    if length(nonzero_rows) > 10
        fprintf('... and %d more\n', length(nonzero_rows) - 10);
    end
else
    fprintf('\nWARNING: No non-zero pixels found!\n');
end

% Plot with different colormaps
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
imagesc(X_mm(1,:), Z_mm(:,1), scene_matrix);
colormap(gray);
colorbar;
axis image;
set(gca, 'YDir', 'normal');
xlabel('X Position (mm)');
ylabel('Z Position (mm)');
title('Gray Colormap');

subplot(1, 3, 2);
imagesc(X_mm(1,:), Z_mm(:,1), scene_matrix);
colormap(jet);
colorbar;
axis image;
set(gca, 'YDir', 'normal');
xlabel('X Position (mm)');
ylabel('Z Position (mm)');
title('Jet Colormap');

subplot(1, 3, 3);
imagesc(X_mm(1,:), Z_mm(:,1), scene_matrix);
colormap(hot);
colorbar;
axis image;
set(gca, 'YDir', 'normal');
xlabel('X Position (mm)');
ylabel('Z Position (mm)');
title('Hot Colormap');

sgtitle(sprintf('Scene Matrix Visualization (%.1f%% non-zero pixels)', nnz(scene_matrix)/numel(scene_matrix)*100));

% Save debug plot
saveas(gcf, 'debug_scene_analysis.png');
close(gcf);

fprintf('\nDebug plot saved as: debug_scene_analysis.png\n'); 