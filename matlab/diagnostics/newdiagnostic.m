% =========================================================================
% Definitive Data Diagnostic Script v2.0 (Downsample & Filter)
%
% Description:
% This version fixes the filtering issue by first downsampling the raw
% data to a more reasonable rate, which allows for stable filter design.
% It also correctly applies the Butterworth filter as requested.
% =========================================================================

clearvars;
clc;
close all;

%% --- Get User Input for File and Settings ---
fprintf('--- Definitive Data Diagnostic Tool ---\n');
MEASURED_DATA_H5_FILE = input('Enter the path to the HDF5 data file: ', 's');
if ~exist(MEASURED_DATA_H5_FILE, 'file')
    error('HDF5 data file not found: %s', MEASURED_DATA_H5_FILE);
end

fprintf('\n--- Please provide the settings used for this acquisition ---\n');
timebase = input('Enter the timebase integer used (e.g., 4): ');

% --- Calculate Sampling Frequency (fs) from Timebase ---
% This formula is from the PicoScope 5000A Series Programmer's Guide
if timebase < 3
    sample_interval_ns = 2^timebase;
else
    sample_interval_ns = (2^(timebase - 3)) * 8;
end
fs = 1 / (sample_interval_ns * 1e-9);
fprintf('Based on timebase %d, the calculated Sampling Frequency (fs) is: %.2f MHz\n', timebase, fs/1e6);

%% Load and Process Experimental Data
fprintf('\n--- Loading and Processing Experimental Data ---\n');
% Inspect the HDF5 file to find the correct dataset name
info = h5info(MEASURED_DATA_H5_FILE);
dataset_name = info.Datasets(1).Name;
fprintf('Found and using dataset: ''%s''\n', dataset_name);
echo_data_from_h5 = h5read(MEASURED_DATA_H5_FILE, ['/' dataset_name]);
% Transpose the matrix to fix the row/column order from Python
echo_data_matrix = echo_data_from_h5';
[R_acquisitions, total_samples_per_acq] = size(echo_data_matrix);
% Assume 10% pre-trigger samples, as set in the Python script
pre_trigger_samples = round(total_samples_per_acq * 0.1);
fprintf('Loaded and transposed data. Found %d acquisitions of %d samples each.\n', R_acquisitions, total_samples_per_acq);

% --- Signal Processing ---
fprintf('Averaging all acquisitions...\n');
u_avg = mean(echo_data_matrix, 1);
u_avg_dc_corrected = u_avg - mean(u_avg);

% *** NEW: Downsample the signal before filtering ***
downsample_factor = 25;
fs_new = fs / downsample_factor;
u_avg_downsampled = decimate(u_avg_dc_corrected, downsample_factor);
fprintf('Downsampled signal from %.2f MHz to %.2f MHz.\n', fs/1e6, fs_new/1e6);

% *** NEW: Design and apply Butterworth filter to the downsampled signal ***
cutoff_freq = 80e3; 
nyquist_freq = fs_new / 2;
filter_order = 4;
[b, a] = butter(filter_order, cutoff_freq / nyquist_freq, 'low');
u_avg_filtered = filtfilt(b, a, u_avg_downsampled);
fprintf('Applied %d-order Butterworth filter with %.0f kHz cutoff.\n', filter_order, cutoff_freq/1e3);

fprintf('Calculating signal envelope of the filtered signal...\n');
signal_envelope = abs(hilbert(u_avg_filtered));

% Create time axes for both raw and filtered data
time_axis_us_raw = ( double(0:total_samples_per_acq-1) - double(pre_trigger_samples) ) / fs * 1e6;
pre_trigger_downsampled = round(pre_trigger_samples / downsample_factor);
time_axis_us_filt = ( double(0:length(u_avg_filtered)-1) - double(pre_trigger_downsampled) ) / fs_new * 1e6;


%% --- Plotting and Analysis ---
fprintf('Plotting results...\n');
figure(1);
clf;
set(gcf, 'Position', [50, 50, 1400, 800], 'Color', 'w');

% --- 1. Plot the FULL Averaged Signal (Raw vs. Filtered) ---
subplot(2, 2, 1);
plot(time_axis_us_raw, u_avg_dc_corrected, 'b', 'DisplayName', 'Raw Averaged Signal');
hold on;
plot(time_axis_us_filt, u_avg_filtered, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered Signal');
hold off;
grid on; axis tight;
xlabel('Time (\mus)');
ylabel('Raw ADC Value');
title('Full Signal: Raw vs. Filtered');
legend;

% --- 2. Plot the FULL Filtered Signal Envelope ---
subplot(2, 2, 2);
plot(time_axis_us_filt, signal_envelope, 'r', 'LineWidth', 1.5);
grid on; axis tight;
xlabel('Time (\mus)');
ylabel('Envelope Amplitude');
title('Full Envelope of Filtered Signal');

% --- 3. Zoom in on the "Main Bang" ---
subplot(2, 2, 3);
plot(time_axis_us_filt, signal_envelope, 'r', 'LineWidth', 1.5);
grid on;
xlim([-10, 50]); % Zoom in on the first 50us
xlabel('Time (\mus)');
ylabel('Envelope Amplitude');
title('Zoom-in: Main Bang & Near-Field');

% --- 4. Zoom in on the Expected Echo Region ---
subplot(2, 2, 4);
plot(time_axis_us_filt, signal_envelope, 'r', 'LineWidth', 1.5);
grid on;
xlim([300, 500]); % Zoom in on the 300-500us region
xlabel('Time (\mus)');
ylabel('Envelope Amplitude');
title('Zoom-in: Expected Echo Region');

sgtitle('Definitive Data Diagnostic (with Filtering)');

%% --- Final Analysis ---
fprintf('\n--- Overall Analysis ---\n');
[~, max_env_idx] = max(signal_envelope);
peak_time_us = time_axis_us_filt(max_env_idx);
fprintf('The absolute peak of the signal envelope occurs at %.2f us.\n', peak_time_us);
fprintf('Examine the plots to identify the main bang and any subsequent echo peaks.\n');
