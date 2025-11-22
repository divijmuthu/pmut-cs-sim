% =========================================================================
% Diagnostic Script v5.1 (Path Fix)
%
% Description:
% This version fixes a bug in the file path generation that caused the
% script to fail when finding the corresponding .csv files. It now
% correctly parses the timestamp from the selected .h5 file.
% =========================================================================
clearvars;
clc;
close all;

%% --- 1. Configuration ---
fprintf('--- 1. Configuring Diagnostic ---\n');

% --- Input Files ---
% Select the H5 file using the dialog. The script will find the rest.
H5_FILE_PATH = ''; % Leave empty to open file dialog, or hardcode a path.

% --- Analysis Parameters ---
ACQS_TO_COMPARE = [1, 8, 19, 27]; % Which acquisitions to compare for incoherence

% --- Filter Parameters ---
FC_CUTOFF = 80e3; % Lowpass filter cutoff frequency in Hz

%% --- 2. Load and Pre-process Data ---
fprintf('\n--- 2. Loading and Pre-processing Data ---\n');

% --- Select H5 file if not hardcoded ---
if ~exist(H5_FILE_PATH, 'file') || isempty(H5_FILE_PATH)
    [h5_filename, data_dir] = uigetfile('*.h5', 'Select the HDF5 Data File');
    if isequal(h5_filename,0)
        disp('User canceled file selection. Exiting.');
        return;
    end
    MEASURED_DATA_H5_FILE = fullfile(data_dir, h5_filename);
else
    MEASURED_DATA_H5_FILE = H5_FILE_PATH;
    [data_dir, ~, ~] = fileparts(MEASURED_DATA_H5_FILE);
end

% *** FIX: Correctly parse the timestamp from the filename ***
[~, base_name_no_ext, ~] = fileparts(MEASURED_DATA_H5_FILE);
timestamp = extractAfter(base_name_no_ext, 'data_');
if isempty(timestamp)
    error('Could not extract a valid timestamp from the filename. Expected format: "...data_YYYY-MM-DD_HH-MM-SS.h5"');
end

DELAY_PROFILE_CSV_FILE = fullfile(data_dir, ['delays_', timestamp, '.csv']);
FREQ_PROFILE_CSV_FILE = fullfile(data_dir, ['frequencies_', timestamp, '.csv']);

if ~exist(DELAY_PROFILE_CSV_FILE, 'file'), error('Could not find matching delays file:\n%s', DELAY_PROFILE_CSV_FILE); end
if ~exist(FREQ_PROFILE_CSV_FILE, 'file'), error('Could not find matching frequencies file:\n%s', FREQ_PROFILE_CSV_FILE); end

fprintf('Found and loaded corresponding files for timestamp: %s\n', timestamp);
delay_profiles_us = readmatrix(DELAY_PROFILE_CSV_FILE);
frequencies_hz = readmatrix(FREQ_PROFILE_CSV_FILE);
[R_acquisitions, ~] = size(delay_profiles_us);

% --- Load H5 Data and Convert to Voltage ---
info = h5info(MEASURED_DATA_H5_FILE);
dataset_name = info.Datasets(1).Name;
echo_data_matrix_raw = h5read(MEASURED_DATA_H5_FILE, ['/' dataset_name])';
try
    adc_max = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'max_adc');
    adc_min = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'min_adc');
    v_max = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'max_volts');
    v_min = h5readatt(MEASURED_DATA_H5_FILE, ['/' dataset_name], 'min_volts');
    m = (v_max - v_min) / (adc_max - adc_min);
    c_offset = v_max - m * adc_max;
    echo_data_matrix_volts = double(echo_data_matrix_raw) * m + c_offset;
catch
    warning('Could not read voltage attributes. Using raw data as voltage.');
    echo_data_matrix_volts = double(echo_data_matrix_raw);
end
clear echo_data_matrix_raw;

% --- Calculate True Sampling Rate ---
[~, total_samples_raw] = size(echo_data_matrix_volts);
total_duration_s = 10e-3;
fs_raw = total_samples_raw / total_duration_s;
fprintf('Original sampling rate: %.2f MHz\n', fs_raw/1e6);

% --- Lowpass Filter ---
[b_filter, a_filter] = butter(4, FC_CUTOFF / (fs_raw/2), 'low');
echo_data_matrix_filtered = zeros(size(echo_data_matrix_volts));
for i = 1:R_acquisitions
    echo_data_matrix_filtered(i, :) = filtfilt(b_filter, a_filter, echo_data_matrix_volts(i, :));
end
clear echo_data_matrix_volts;

%% --- 3. Analyze and Visualize Correlation ---
fprintf('\n--- 3. Analyzing Correlation ---\n');

figure('Position', [100, 100, 1400, 700]);
time_axis_ms = (0:total_samples_raw-1)/fs_raw * 1000;

% Ensure we only try to plot acquisitions that exist
acqs_to_plot_safe = ACQS_TO_COMPARE(ACQS_TO_COMPARE <= R_acquisitions);
if isempty(acqs_to_plot_safe) && R_acquisitions > 0
    acqs_to_plot_safe = 1; % Default to plotting the first one
end
fprintf('Comparing acquisitions: %s\n', mat2str(acqs_to_plot_safe));

% --- Plot the selected acquisitions ---
hold on;
for i = 1:length(acqs_to_plot_safe)
    acq_idx = acqs_to_plot_safe(i);
    freq_khz = frequencies_hz(acq_idx) / 1000;
    plot(time_axis_ms, echo_data_matrix_filtered(acq_idx, :), ...
        'DisplayName', sprintf('Acq #%d (%.2f kHz)', acq_idx, freq_khz));
end
hold off;
grid on;
axis tight;
legend('Location', 'northeast');
title('Visual Comparison of Different Acquisitions');
xlabel('Time (ms)');
ylabel('Filtered Amplitude (V)');
xlim([0 2]); % Zoom in on the first 2ms to see the echo clearly

% --- Calculate and Display Correlation Coefficients ---
if length(acqs_to_plot_safe) > 1
    ref_acq_idx = acqs_to_plot_safe(1);
    ref_signal = echo_data_matrix_filtered(ref_acq_idx, :);
    
    corr_string = {'Correlation Coefficients:', sprintf('(vs. Acq #%d)', ref_acq_idx), ''};
    fprintf('\nCorrelation Coefficients (vs. Acq #%d):\n', ref_acq_idx);
    
    for i = 2:length(acqs_to_plot_safe)
        curr_acq_idx = acqs_to_plot_safe(i);
        curr_signal = echo_data_matrix_filtered(curr_acq_idx, :);
        
        corr_matrix = corrcoef(ref_signal, curr_signal);
        corr_val = corr_matrix(1, 2);
        
        fprintf('  - Acq #%d: %.4f\n', curr_acq_idx, corr_val);
        corr_string{end+1} = sprintf('Acq #%d: %.4f', curr_acq_idx, corr_val);
    end
    
    annotation('textbox', [0.75, 0.6, 0.15, 0.2], 'String', corr_string, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
end

disp('Diagnostic script finished.');
