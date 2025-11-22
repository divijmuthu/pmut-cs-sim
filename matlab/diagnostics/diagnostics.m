% =========================================================================
% Comprehensive Data Diagnostic Script v4.1 (Lowpass Filtered)
%
% Description:
% This version uses a 4th-order Butterworth lowpass filter with an 80 kHz
% cutoff to align with previous analysis. The filter is applied to both
% measured and simulated data for an accurate comparison.
% =========================================================================
clearvars;
clc;
close all;

% --- Configuration ---
% Point to the data files you want to analyze
DELAY_PROFILE_CSV_FILE = 'delays_2025-07-07_17-03-00.csv';
MEASURED_DATA_H5_FILE = 'data_2025-07-07_17-03-00.h5';

% --- Analysis Parameters ---
ACQS_TO_ANALYZE = [1, 18, 42]; % Which acquisitions to compare for incoherence
ACQUISITION_FOR_SIM_COMPARISON = 1; % Which acquisition to use for the detailed comparison

% --- Initialize Field II ---
field_init(-1);

% --- Core Physical and Simulation Constants ---
c = 1540;                   % Speed of Sound [m/s]
fs = 125e6;                 % The TRUE sampling rate of the collected data
set_field('fs', fs);
set_field('c', c);

fprintf('--- Comprehensive Diagnostic Script v4.1 (Lowpass Filtered) ---\n');
fprintf('Analyzing data at the original %.1f MHz sampling rate.\n\n', fs/1e6);

% --- Geometry and Simulation Parameters ---
pMUT_width_mm = 20; pMUT_spacing_mm = 25; kerf_mm = 0.1;
excitation_amplitude = 1; % Field II works with normalized pressure, so we use 1
% --- Convert mm to meters ---
pMUT_width = pMUT_width_mm/1000; pMUT_height = pMUT_width; kerf = kerf_mm/1000;
d_spacing = pMUT_spacing_mm / 1000;

% --- Define a SINGLE target point for simulation ---
target_pos = [0, 0, 0.15]; % A single point target at 15cm depth

% --- Define Separate Tx and Rx Apertures using Grid Mapping ---
fprintf('Defining separate Tx and Rx apertures...\n');
radius = d_spacing;
tx_positions = [radius, 0, 0; radius*cos(2*pi/3), radius*sin(2*pi/3), 0; radius*cos(4*pi/3), radius*sin(4*pi/3), 0];
rx_position = [0, 0, 0];
num_x_grid = 15; num_y_grid = 15;
physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
element_no_grid_map = 0;
center_offset_x = (num_x_grid - 1) / 2 * (pMUT_width + kerf);
center_offset_y = (num_y_grid - 1) / 2 * (pMUT_height + kerf);
for iy = 1:num_y_grid, y_pos_el = (iy - 1) * (pMUT_height + kerf) - center_offset_y;
for ix = 1:num_x_grid, x_pos_el = (ix - 1) * (pMUT_width + kerf) - center_offset_x;
element_no_grid_map = element_no_grid_map + 1;
physical_element_centers(element_no_grid_map, :) = [x_pos_el, y_pos_el, 0]; end, end
tx_indices_linear = zeros(size(tx_positions,1), 1);
for i = 1:size(tx_positions,1), [~, tx_indices_linear(i)] = min(sum((physical_element_centers - tx_positions(i,:)).^2, 2)); end
tx_indices_linear = unique(tx_indices_linear); num_active_tx = length(tx_indices_linear);
tx_enabled_matrix = zeros(num_y_grid, num_x_grid); tx_enabled_matrix(tx_indices_linear) = 1;
Tx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, tx_enabled_matrix, 1, 1, [0 0 0]);
[~, rx_index_linear] = min(sum((physical_element_centers - rx_position).^2, 2));
rx_enabled_matrix = zeros(num_y_grid, num_x_grid); rx_enabled_matrix(rx_index_linear) = 1;
Rx_Aperture = xdc_2d_array(num_x_grid, num_y_grid, pMUT_width, pMUT_height, kerf, kerf, rx_enabled_matrix, 1, 1, [0 0 0]);
num_active = num_active_tx;
fprintf('Successfully mapped %d Tx and 1 Rx element.\n', num_active);

% --- Define Impulse Response ---
f_start_chirp = 10e3; f_end_chirp = 200e3; burst_duration = 0.02e-3;
t_burst_vec = 0 : 1/fs : burst_duration;
synth_burst_base = chirp(t_burst_vec, f_start_chirp, t_burst_vec(end), f_end_chirp, 'linear');
synth_burst_windowed = synth_burst_base .* tukeywin(length(t_burst_vec), 0.25)';
impulse_response_waveform = synth_burst_windowed * excitation_amplitude;
xdc_impulse(Tx_Aperture, impulse_response_waveform); xdc_excitation(Tx_Aperture, 1);

%% Load and Process Experimental Data
fprintf('\n--- Loading and Processing Experimental Data ---\n');
if ~exist(DELAY_PROFILE_CSV_FILE, 'file'), error('Delay CSV not found!'); end
delay_profiles_us = readmatrix(DELAY_PROFILE_CSV_FILE);

if ~exist(MEASURED_DATA_H5_FILE, 'file'), error('HDF5 data file not found!'); end
info = h5info(MEASURED_DATA_H5_FILE);
dataset_name = info.Datasets(1).Name;

% --- ADC to Voltage Conversion ---
echo_data_matrix_raw = h5read(MEASURED_DATA_H5_FILE, ['/' dataset_name])';
try
    adc_max = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'max_adc');
    adc_min = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'min_adc');
    v_max = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'max_volts');
    v_min = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'min_volts');
    m = (v_max - v_min) / (adc_max - adc_min);
    c_offset = v_max - m * adc_max;
    echo_data_matrix_volts = double(echo_data_matrix_raw) * m + c_offset;
    fprintf('Successfully converted raw ADC values to Voltage.\n');
catch ME
    warning('Could not find voltage conversion attributes. Using raw ADC values.');
    fprintf('Error: %s\n', ME.message);
    echo_data_matrix_volts = double(echo_data_matrix_raw);
end

[R_acquisitions, total_samples] = size(echo_data_matrix_volts);
pre_trigger_samples = round(total_samples * 0.1);
fprintf('Loaded data for %d acquisitions of %d samples each.\n', R_acquisitions, total_samples);

%% --- NEW: Apply Butterworth Lowpass Filter ---
fprintf('\n--- Applying 4th-Order Butterworth Lowpass Filter ---\n');
N_order = 4;                % Filter order
fc_cutoff = 80e3;           % Cutoff frequency (80 kHz)
Wn = fc_cutoff / (fs/2);    % Normalize frequency to Nyquist
[b, a] = butter(N_order, Wn, 'low'); % Get filter coefficients for LOWPASS

fprintf('Filter Cutoff: %.1f kHz\n', fc_cutoff/1e3);

% Apply the filter to each acquisition (each row)
echo_data_matrix_filtered = zeros(size(echo_data_matrix_volts));
for i = 1:R_acquisitions
    % Use filtfilt for zero-phase distortion
    echo_data_matrix_filtered(i, :) = filtfilt(b, a, echo_data_matrix_volts(i, :));
end
fprintf('Filtering complete.\n');


%% --- Diagnostic 1: Visual Comparison & Correlation of FILTERED Data ---
fprintf('\n--- Diagnostic 1: Verifying Incoherence (on Filtered Data) ---\n');
figure(1); clf; set(gcf, 'Position', [100, 100, 1400, 600], 'Color', 'w');
time_axis_us = ( (0:total_samples-1) - pre_trigger_samples) / fs * 1e6;

hold on;
acqs_to_plot_safe = ACQS_TO_ANALYZE(ACQS_TO_ANALYZE <= R_acquisitions);
if isempty(acqs_to_plot_safe) && R_acquisitions > 0
    acqs_to_plot_safe = 1;
end
fprintf('Displaying acquisitions: %s\n', mat2str(acqs_to_plot_safe));

for i = 1:length(acqs_to_plot_safe)
    acq_idx = acqs_to_plot_safe(i);
    % *** USE FILTERED DATA ***
    plot(time_axis_us, echo_data_matrix_filtered(acq_idx, :), 'DisplayName', sprintf('Acq #%d', acq_idx));
end
hold off; grid on; axis tight; legend('Location', 'northeast');
title('Visual Comparison of Different Acquisitions (Filtered)');
xlabel('Time (\mus)');
ylabel('Voltage (V)');

% Correlation Check on FILTERED data
corr_coeffs = [];
if length(acqs_to_plot_safe) > 1
    ref_signal = echo_data_matrix_filtered(acqs_to_plot_safe(1), :);
    corr_coeffs = zeros(length(acqs_to_plot_safe), 1);
    for i = 1:length(acqs_to_plot_safe)
        acq_idx = acqs_to_plot_safe(i);
        corr_matrix = corrcoef(ref_signal, echo_data_matrix_filtered(acq_idx, :));
        corr_coeffs(i) = corr_matrix(1, 2);
        fprintf('Correlation between Acq #%d and Acq #%d: %.4f\n', acqs_to_plot_safe(1), acq_idx, corr_coeffs(i));
    end
    
    corr_string = {'Correlation Coefficients:', ''};
    for i = 2:length(acqs_to_plot_safe)
        corr_string{end+1} = sprintf('Acq #%d vs #%d: %.4f', acqs_to_plot_safe(1), acqs_to_plot_safe(i), corr_coeffs(i));
    end
    annotation('textbox', [0.7, 0.5, 0.2, 0.2], 'String', corr_string, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');
end

%% --- Diagnostic 2: Measured vs. Simulated (FILTERED) ---
fprintf('\n--- Diagnostic 2: Comparing Measured vs. Ideal Simulated Echo (Filtered) ---\n');
% Simulate one "perfect" echo
delay_vector_s = delay_profiles_us(ACQUISITION_FOR_SIM_COMPARISON, :) / 1e6;
xdc_focus_times(Tx_Aperture, 0, delay_vector_s);
[hhp_sim, start_time_sim] = calc_hhp(Tx_Aperture, Rx_Aperture, target_pos);

t_target_axis = time_axis_us / 1e6;
t_sim_axis = start_time_sim + (0:(length(hhp_sim)-1))/fs;
u_simulated_raw = interp1(t_sim_axis, hhp_sim, t_target_axis, 'linear', 0);

% *** USE FILTERED DATA ***
u_measured_single_filtered = echo_data_matrix_filtered(ACQUISITION_FOR_SIM_COMPARISON, :);
u_measured_single_filtered = u_measured_single_filtered - mean(u_measured_single_filtered);

% *** APPLY SAME FILTER TO SIMULATED DATA FOR FAIR COMPARISON ***
u_simulated_filtered = filtfilt(b, a, u_simulated_raw);

% --- Automatic Amplitude Scaling (on filtered signals) ---
rms_measured = rms(u_measured_single_filtered);
rms_simulated = rms(u_simulated_filtered);
scaling_factor = rms_measured / rms_simulated;
u_simulated_scaled_filtered = u_simulated_filtered * scaling_factor;

% --- Plot Comparison ---
figure(2); clf; set(gcf, 'Position', [150, 150, 1400, 800], 'Color', 'w');

% Envelopes of filtered signals
env_measured = abs(hilbert(u_measured_single_filtered));
env_simulated_scaled = abs(hilbert(u_simulated_scaled_filtered));

subplot(2,1,1);
plot(time_axis_us, u_measured_single_filtered, 'b', 'DisplayName', 'Measured (Filtered)');
hold on;
plot(time_axis_us, u_simulated_scaled_filtered, 'r--', 'DisplayName', 'Simulated (Filtered & Scaled)');
hold off; grid on; axis tight; legend;
title('Waveform Comparison');
ylabel('Amplitude (V)');

subplot(2,1,2);
plot(time_axis_us, env_measured, 'b', 'LineWidth', 1.5, 'DisplayName', 'Measured Envelope');
hold on;
plot(time_axis_us, env_simulated_scaled, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Simulated Envelope (Scaled)');
hold off; grid on; axis tight; legend;
title('Envelope Comparison');
xlabel('Time (\mus)');
ylabel('Envelope Amplitude (V)');

sgtitle(sprintf('Diagnostic Comparison for Acquisition #%d (Filtered)', ACQUISITION_FOR_SIM_COMPARISON));

% --- Frequency Content (FFT) Comparison (of filtered signals) ---
figure(3); clf; set(gcf, 'Position', [200, 200, 1200, 600], 'Color', 'w');
NFFT = 2^nextpow2(total_samples);
Y_measured = fft(u_measured_single_filtered, NFFT)/total_samples;
Y_simulated = fft(u_simulated_scaled_filtered, NFFT)/total_samples;
f_axis = fs/2*linspace(0,1,NFFT/2+1);

plot(f_axis/1000, 2*abs(Y_measured(1:NFFT/2+1)), 'b', 'DisplayName', 'Measured Signal (Filtered)');
hold on;
plot(f_axis/1000, 2*abs(Y_simulated(1:NFFT/2+1)), 'r--', 'DisplayName', 'Simulated Signal (Filtered & Scaled)');
hold off;
title('Frequency Spectrum Comparison');
xlabel('Frequency (kHz)');
ylabel('Amplitude (V)');
grid on; xlim([0 250]); legend;

%% End Field II Simulation
if exist('field_end', 'file') == 2, field_end; disp('Field II ended.'); end
