% =========================================================================
% Interactive Echo Finder & Diagnostic Script v10 (Smart Dataset Selection)
%
% Description:
% This version is now robust to different HDF5 file structures. It inspects
% the file, lists the available datasets, and asks the user to choose
% which one to analyze.
% =========================================================================

clearvars;
clc;
close all;

% --- Configuration ---
MEASURED_DATA_H5_FILE = 'data_2025-07-07_17-03-00.h5';
PLOT_CHUNK_SIZE_US = 100; % Display 100us at a time

%% --- Get User Input for Settings ---
fprintf('--- Interactive Data Verification ---\n');
fprintf('Using data file: %s\n', MEASURED_DATA_H5_FILE);
if ~exist(MEASURED_DATA_H5_FILE, 'file')
    error('HDF5 data file not found: %s', MEASURED_DATA_H5_FILE);
end

% *** UPGRADE: Inspect the HDF5 file and ask user to choose the dataset ***
try
    info = h5info(MEASURED_DATA_H5_FILE);
    fprintf('\nAvailable datasets in this file:\n');
    for i = 1:length(info.Datasets)
        fprintf('  - %s\n', info.Datasets(i).Name);
    end
    dataset_name = input('\nPlease type the name of the echo data dataset to analyze: ', 's');
catch ME
    error('Could not inspect HDF5 file. Error: %s', ME.message);
end


fprintf('\n--- Please provide the settings used for this acquisition ---\n');
timebase = input('Enter the timebase integer used (e.g., 6): ');

% Calculate sampling frequency based on the timebase
if timebase < 3
    fs = 1 / (2^timebase * 1e-9);
else
    fs = 125e6 / (2^(timebase - 3));
end
fprintf('Calculated Sampling Frequency (fs): %.2f MHz\n', fs/1e6);

%% Load Experimental Data
fprintf('\n--- Loading Experimental Data ---\n');
try
    echo_data_from_h5 = h5read(MEASURED_DATA_H5_FILE, ['/' dataset_name]);
catch ME
    error('Could not read dataset "%s". Please check the name and try again. Error: %s', dataset_name, ME.message);
end

echo_data_matrix = echo_data_from_h5'; % Transpose to get (Acqs x Samples)
[R_acquisitions, total_samples_per_acq] = size(echo_data_matrix);
pre_trigger_samples = round(total_samples_per_acq * 0.1); % Assume 10% pre-trigger
fprintf('Loaded data for %d acquisitions of %d samples each.\n', R_acquisitions, total_samples_per_acq);

%% --- Signal Processing to Find Echo ---
fprintf('Averaging all acquisitions...\n');
u_avg = mean(echo_data_matrix, 1);
u_avg_dc_corrected = u_avg - mean(u_avg);
fprintf('Calculating signal envelope...\n');
signal_envelope = abs(hilbert(u_avg_dc_corrected));
time_axis_us = ( double(0:total_samples_per_acq-1) - double(pre_trigger_samples) ) / fs * 1e6;

%% --- Plotting and Analysis ---
fprintf('Plotting results in segments...\n');
figure(1);
clf;
set(gcf, 'Position', [100, 100, 1400, 700], 'Color', 'w');

num_chunks = ceil(time_axis_us(end) / PLOT_CHUNK_SIZE_US);
if isempty(num_chunks) || num_chunks == 0, num_chunks = 1; end
t_start = time_axis_us(1);

for i = 1:num_chunks
    t_end = t_start + PLOT_CHUNK_SIZE_US;
    
    idx_start = find(time_axis_us >= t_start, 1, 'first');
    idx_end = find(time_axis_us <= t_end, 1, 'last');
    
    if isempty(idx_start), idx_start = 1; end
    if isempty(idx_end), idx_end = length(time_axis_us); end
    
    % Plot the averaged signal for this chunk
    subplot(2, 1, 1);
    plot(time_axis_us(idx_start:idx_end), u_avg_dc_corrected(idx_start:idx_end));
    grid on;
    xlim([t_start, t_end]);
    xlabel('Time (\mus)');
    ylabel('Raw ADC Value');
    title(sprintf('Averaged Signal (Segment %d/%d)', i, num_chunks));

    % Plot the signal envelope for this chunk
    subplot(2, 1, 2);
    plot(time_axis_us(idx_start:idx_end), signal_envelope(idx_start:idx_end), 'r', 'LineWidth', 1.5);
    grid on;
    xlim([t_start, t_end]);
    xlabel('Time (\mus)');
    ylabel('Envelope Amplitude');
    title('Signal Envelope (Look for a Peak)');
    
    sgtitle(sprintf('Diagnostic Scan: %.1f us to %.1f us', t_start, t_end));
    
    fprintf('Displaying time segment %.1f us to %.1f us. Press any key to continue...\n', t_start, t_end);
    pause; % Wait for user to press a key
    
    t_start = t_end;
    if t_start > time_axis_us(end)
        break;
    end
end

fprintf('\n--- Overall Analysis ---\n');
[~, max_env_idx] = max(signal_envelope);
peak_time_us = time_axis_us(max_env_idx);
fprintf('The absolute peak of the signal envelope occurs at %.2f us.\n', peak_time_us);
fprintf('If you saw a distinct peak in one of the segments, that is your echo.\n');

