% =========================================================================
% Echo Finder & Diagnostic Script v7 (Transpose Fix)
%
% Description:
% This version fixes a critical bug by transposing the data matrix after
% loading it from HDF5. This corrects the row-major/column-major mismatch
% between Python and MATLAB.
% =========================================================================

clearvars; 
clc;
close all;

% --- Configuration: Point to the data file you want to inspect ---
MEASURED_DATA_H5_FILE = 'data_2025-07-03_17-13-37.h5';

% --- Core Constants ---
fs = 125e6; % Sampling Frequency [Hz]
PLOT_CHUNK_SIZE_US = 100; % Display 100us at a time

%% Load Experimental Data
fprintf('\n--- Loading Experimental Data to Find Echo ---\n');
if ~exist(MEASURED_DATA_H5_FILE, 'file')
    error('HDF5 data file not found: %s', MEASURED_DATA_H5_FILE);
end

% Load the data and metadata
echo_data_from_h5 = h5read(MEASURED_DATA_H5_FILE, '/echo_data');

% *** FIX APPLIED HERE: Transpose the matrix to fix row/column order ***
echo_data_matrix = echo_data_from_h5';

[R_acquisitions, total_samples_per_acq] = size(echo_data_matrix);
try
    pre_trigger_samples = h5readatt(MEASURED_DATA_H5_FILE, '/', 'pre_trigger_samples');
catch
    pre_trigger_samples = round(total_samples_per_acq * 0.1);
end
fprintf('Loaded and transposed data. Found %d acquisitions of %d samples each.\n', R_acquisitions, total_samples_per_acq);

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
    
    % Find indices for the current chunk
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
[max_env_val, max_env_idx] = max(signal_envelope);
peak_time_us = time_axis_us(max_env_idx);
fprintf('The absolute peak of the signal envelope occurs at %.2f us.\n', peak_time_us);
fprintf('If you saw a distinct peak in one of the segments, that is your echo.\n');

