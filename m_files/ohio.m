% Diagnostic_Signal_Test_v1.m
%
% GOAL: A radically simplified script to test ONE thing:
% Can we generate a strong H matrix by defining the chirp as the
% transducer's fundamental impulse response, bypassing ele_waveform?
%
% This is NOT a CS script. It's a test to find the source of the
% low signal amplitude.
%
clearvars; clc; close all;

% --- Initialize Field II ---
field_init(-1);

% --- Use the Corrected Physical Scaling ---
c = 1540;
fs = 2e6;
fc_nominal = 1.0e5; % <<<< Using the corrected 100 kHz
lambda = c/fc_nominal;
fprintf('Running Diagnostic with fc_nominal = %g kHz\n', fc_nominal/1e3);
set_field('fs', fs);
set_field('c', c);

% --- Define Synthetic Chirp ---
% This will now be our impulse response.
f_start_chirp = 10e3;
f_end_chirp = 200e3;
chirp_duration = 0.2e-3;
t_chirp_vec = 0 : 1/fs : chirp_duration;
synth_chirp_base = chirp(t_chirp_vec, f_start_chirp, t_chirp_vec(end), f_end_chirp, 'linear');
synth_chirp_windowed = synth_chirp_base .* tukeywin(length(t_chirp_vec), 0.25)';

% --- Define Transducer Impulse Response and Excitation ---
% RADICAL CHANGE: The chirp IS the impulse response.
impulseResponse = synth_chirp_windowed;
excitationPulse = [1]; % Trigger with a single delta function.
fprintf('DIAGNOSTIC: Using the chirp as the main impulse response.\n');

% --- Transducer and Grid Definition (Same as v3/v4) ---
pMUT_width = 1.5*lambda;
pMUT_height = pMUT_width;
kerf = 0.5*lambda;
d_spacing = 5*lambda;
pos1 = [d_spacing,0,0];
pos2 = [d_spacing*cos(2*pi/3), d_spacing*sin(2*pi/3),0];
pos3 = [d_spacing*cos(4*pi/3), d_spacing*sin(4*pi/3),0];
desired_positions = [pos1; pos2; pos3];
num_active_intended = size(desired_positions, 1);
array_span_x = 2.5*d_spacing;
array_span_y = 2.5*d_spacing;
element_width_grid = pMUT_width;
element_height_grid = pMUT_height;
num_x_grid = ceil(array_span_x / (element_width_grid + kerf));
num_y_grid = ceil(array_span_y / (element_height_grid + kerf));
if mod(num_x_grid, 2) == 0; num_x_grid = num_x_grid + 1; end
if mod(num_y_grid, 2) == 0; num_y_grid = num_y_grid + 1; end

physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
element_no_grid_map = 0;
center_offset_x = (num_x_grid - 1)/2 * (element_width_grid + kerf);
center_offset_y = (num_y_grid - 1)/2 * (element_height_grid + kerf);
for iy = 1:num_y_grid
    y_pos_el = (iy-1)*(element_height_grid + kerf) - center_offset_y;
    for ix = 1:num_x_grid
        x_pos_el = (ix-1)*(element_width_grid + kerf) - center_offset_x;
        element_no_grid_map = element_no_grid_map + 1;
        physical_element_centers(element_no_grid_map, :) = [x_pos_el, y_pos_el, 0];
    end
end
active_indices_linear = zeros(num_active_intended, 1);
for i = 1:num_active_intended
    distances = sqrt(sum((physical_element_centers - desired_positions(i,:)).^2, 2));
    [~, min_idx] = min(distances);
    active_indices_linear(i) = min_idx;
end
active_indices_linear = unique(active_indices_linear);
num_active = length(active_indices_linear);
enabled_matrix = zeros(num_y_grid, num_x_grid);
[row_indices, col_indices] = ind2sub([num_y_grid, num_x_grid], active_indices_linear);
for i = 1:num_active
    enabled_matrix(row_indices(i), col_indices(i)) = 1;
end

pMUT_Aperture = xdc_2d_array(num_x_grid, num_y_grid, element_width_grid, element_height_grid, kerf, kerf, enabled_matrix, 1, 1, [0 0 100]);
xdc_impulse(pMUT_Aperture, impulseResponse);
xdc_excitation(pMUT_Aperture, excitationPulse);

% --- Imaging Grid Definition ---
grid_width = 20*lambda;
grid_depth_start = 10*lambda;
grid_depth_end = 40*lambda;
grid_step = lambda/2;
x_coords_img = -grid_width/2 : grid_step : grid_width/2;
z_coords_img = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
Y_mesh = zeros(size(X_mesh));
N_pixels = numel(X_mesh);
hydrophone_positions_img = [X_mesh(:), Y_mesh(:), Z_mesh(:)];

%% --- Calculate H Matrix ---
fprintf('\n--- DIAGNOSTIC: Calculating H for a single, un-coded acquisition ---\n');

% Calculate pulse-echo field. This is our H matrix for one acquisition.
[H, start_time] = calc_hhp(pMUT_Aperture, pMUT_Aperture, hydrophone_positions_img);

if ~isempty(H)
    max_H_val = max(abs(H(:)));
    fprintf('\n********************************************************\n');
    fprintf('*** RESULT: Maximum absolute value in H is: %g   ***\n', max_H_val);
    fprintf('********************************************************\n\n');
else
    fprintf('\n*** ERROR: H matrix is empty!   ***\n\n');
end

% --- Plot a column of H to check its magnitude ---
figure(1); clf;
plot(H(:, round(N_pixels/2)));
title(sprintf('Diagnostic H Matrix Column (Max Val = %g)', max_H_val));
xlabel('Time Sample Index');
ylabel('Amplitude');
grid on;

%% End Field II
if exist('field_end','file') == 2
    field_end;
    disp('Field II ended.');
end