% Simple H Matrix Test
% Test H matrix generation with included functions
clearvars;
clc;
close all;

%% Test parameters
fprintf('--- Simple H Matrix Test ---\n');

% Set up simple test parameters
params = struct();
params.c = 343;                    % Speed of Sound [m/s]
params.pMUT_width_mm = 20;         % pMUT width/height (mm)
params.pMUT_spacing_mm = 20;       % Triangle sides (mm)
params.grid_depth_start_mm = 250;  % Start depth (mm)
params.grid_depth_end_mm = 350;    % End depth (mm)
params.grid_width_mm = 150;        % Grid width (mm)
params.grid_step_mm = 2.0;         % Grid step (mm)

% Create simple test data
delay_profiles = [1, 2, 3; 4, 5, 6; 7, 8, 9];  % 3 acquisitions, 3 delays each
frequencies = [57000, 57000, 57000; 58000, 58000, 58000; 59000, 59000, 59000];  % 3 acquisitions, 3 freqs each

fprintf('Test parameters:\n');
fprintf('  Grid: %.0f-%.0f mm depth, %.0f mm width, %.1f mm step\n', ...
    params.grid_depth_start_mm, params.grid_depth_end_mm, params.grid_width_mm, params.grid_step_mm);
fprintf('  Delay profiles: %s\n', mat2str(size(delay_profiles)));
fprintf('  Frequencies: %s\n', mat2str(size(frequencies)));

%% Generate grid
fprintf('\n--- Generating Grid ---\n');

grid_depth_range_mm = params.grid_depth_end_mm - params.grid_depth_start_mm;
grid_depth_steps = round(grid_depth_range_mm / params.grid_step_mm) + 1;
grid_width_steps = round(params.grid_width_mm / params.grid_step_mm) + 1;

x_grid = linspace(-params.grid_width_mm/2, params.grid_width_mm/2, grid_width_steps);
z_grid = linspace(params.grid_depth_start_mm, params.grid_depth_end_mm, grid_depth_steps);
[X, Z] = meshgrid(x_grid, z_grid);

grid_points = [X(:), Z(:)];
num_grid_points = size(grid_points, 1);

fprintf('Grid generation:\n');
fprintf('  Grid dimensions: %d x %d = %d points\n', grid_depth_steps, grid_width_steps, num_grid_points);
fprintf('  X range: [%.1f, %.1f] mm\n', min(grid_points(:, 1)), max(grid_points(:, 1)));
fprintf('  Z range: [%.1f, %.1f] mm\n', min(grid_points(:, 2)), max(grid_points(:, 2)));

%% Generate pMUT positions
fprintf('\n--- Generating pMUT Positions ---\n');

num_elements = size(delay_profiles, 2);

% Simple triangular array
radius = params.pMUT_spacing_mm;
angles = linspace(0, 2*pi, num_elements + 1);
angles = angles(1:end-1); % Remove duplicate 0/2pi

pMUT_positions = zeros(num_elements, 2);
for i = 1:num_elements
    pMUT_positions(i, :) = [radius * cos(angles(i)), radius * sin(angles(i))];
end

fprintf('pMUT positions:\n');
for i = 1:size(pMUT_positions, 1)
    fprintf('  Element %d: (%.1f, %.1f) mm\n', i, pMUT_positions(i, 1), pMUT_positions(i, 2));
end

%% Test single H row calculation
fprintf('\n--- Testing Single H Row Calculation ---\n');

% Test first acquisition
acq_idx = 1;
delays = delay_profiles(acq_idx, :);
freqs = frequencies(acq_idx, :);

fprintf('Testing acquisition %d:\n', acq_idx);
fprintf('  Delays: [%.1f, %.1f, %.1f] us\n', delays);
fprintf('  Frequencies: [%.0f, %.0f, %.0f] Hz\n', freqs);

% Calculate H row manually
H_row = zeros(1, num_grid_points);

for grid_idx = 1:num_grid_points
    grid_point = grid_points(grid_idx, :);
    total_response = 0;
    
    for elem_idx = 1:num_elements
        pMUT_pos = pMUT_positions(elem_idx, :);
        
        % Calculate distance
        distance = norm(grid_point - pMUT_pos);
        
        % Calculate time delay
        time_delay = distance / (params.c * 1e-3);
        element_delay = delays(elem_idx) * 1e-6;
        time_delay = time_delay + element_delay;
        
        % Calculate response
        freq = freqs(elem_idx);
        wavelength = params.c / freq;
        
        spatial_response = exp(-distance / (wavelength * 2));
        phase_factor = cos(2 * pi * distance / wavelength);
        response = spatial_response * phase_factor * exp(-time_delay * 1e3) * 1e-1;
        
        total_response = total_response + response;
    end
    
    H_row(grid_idx) = total_response;
end

% Normalize the row
if norm(H_row) > 0
    H_row = H_row / norm(H_row);
end

fprintf('H row statistics:\n');
fprintf('  Size: %s\n', mat2str(size(H_row)));
fprintf('  Mean: %.6f\n', mean(H_row));
fprintf('  Std: %.6f\n', std(H_row));
fprintf('  Min: %.6f\n', min(H_row));
fprintf('  Max: %.6f\n', max(H_row));
fprintf('  Non-zero elements: %d/%d (%.2f%%)\n', sum(H_row ~= 0), length(H_row), 100*sum(H_row ~= 0)/length(H_row));
fprintf('  Norm: %.6f\n', norm(H_row));

%% Test individual response calculation
fprintf('\n--- Testing Individual Response Calculation ---\n');

% Test a few grid points
test_grid_points = [1, round(num_grid_points/2), num_grid_points];
for test_idx = 1:length(test_grid_points)
    grid_idx = test_grid_points(test_idx);
    grid_point = grid_points(grid_idx, :);
    
    fprintf('Grid point %d (%.1f, %.1f) mm:\n', grid_idx, grid_point(1), grid_point(2));
    
    total_response = 0;
    for elem_idx = 1:num_elements
        pMUT_pos = pMUT_positions(elem_idx, :);
        
        % Calculate distance
        distance = norm(grid_point - pMUT_pos);
        
        % Calculate time delay
        time_delay = distance / (params.c * 1e-3);
        element_delay = delays(elem_idx) * 1e-6;
        time_delay = time_delay + element_delay;
        
        % Calculate response
        freq = freqs(elem_idx);
        wavelength = params.c / freq;
        
        spatial_response = exp(-distance / (wavelength * 2));
        phase_factor = cos(2 * pi * distance / wavelength);
        response = spatial_response * phase_factor * exp(-time_delay * 1e3) * 1e-1;
        
        fprintf('  Element %d: distance=%.1f mm, delay=%.6f s, response=%.6f\n', ...
            elem_idx, distance, time_delay, response);
        
        total_response = total_response + response;
    end
    
    fprintf('  Total response: %.6f\n', total_response);
end

%% Plot results
figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Plot grid and pMUT positions
subplot(2, 3, 1);
scatter(grid_points(:, 1), grid_points(:, 2), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
scatter(pMUT_positions(:, 1), pMUT_positions(:, 2), 100, 'r', 'filled');
title('Grid and pMUT Positions');
xlabel('X (mm)'); ylabel('Z (mm)');
legend('Grid Points', 'pMUT Elements', 'Location', 'best');
grid on;
axis equal;

% Plot H row
subplot(2, 3, 2);
plot(H_row, 'b-', 'LineWidth', 2);
title('H Row Values');
xlabel('Grid Point Index'); ylabel('Response Value');
grid on;

% Plot H row histogram
subplot(2, 3, 3);
histogram(H_row, 50);
title('H Row Distribution');
xlabel('Response Value'); ylabel('Count');
grid on;

% Plot H row as image
subplot(2, 3, 4);
H_row_image = reshape(H_row, grid_depth_steps, grid_width_steps);
imagesc(H_row_image);
colorbar;
title('H Row as Image');
xlabel('X (pixels)'); ylabel('Z (pixels)');
axis equal;

% Plot response vs distance
subplot(2, 3, 5);
distances = zeros(num_grid_points, 1);
for i = 1:num_grid_points
    distances(i) = norm(grid_points(i, :) - pMUT_positions(1, :));
end
scatter(distances, H_row, 'b.');
title('Response vs Distance to Element 1');
xlabel('Distance (mm)'); ylabel('Response Value');
grid on;

% Plot response vs grid position
subplot(2, 3, 6);
scatter(grid_points(:, 1), grid_points(:, 2), 20, H_row, 'filled');
colorbar;
title('Response vs Grid Position');
xlabel('X (mm)'); ylabel('Z (mm)');
axis equal;

sgtitle('Simple H Matrix Test');
set(gcf, 'Color', 'w');

fprintf('\n--- Test Complete ---\n');
fprintf('Check the plots to see if the H row calculation is working!\n'); 