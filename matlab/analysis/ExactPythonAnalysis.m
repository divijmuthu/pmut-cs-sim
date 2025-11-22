% Exact Python Analysis Replication
% This script exactly replicates the Python SNR calculation to identify discrepancies
clearvars;
clc;
close all;

%% ===== EXACT PYTHON REPLICATION =====
fprintf('--- Exact Python Analysis Replication ---\n');

% Load the same data
h5_file_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/data_2025-07-18_15-08-13_subtracted.h5';
freq_csv_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/frequencies_2025-07-18_15-37-31.csv';

% Load raw ADC data (exactly like Python)
raw_adc_data = h5read(h5_file_path, '/echo_data_raw_adc');
if size(raw_adc_data, 1) > size(raw_adc_data, 2)
    raw_adc_data = raw_adc_data';
end

% Load frequencies
frequencies_hz = readmatrix(freq_csv_path);

fprintf('Raw ADC data shape: %s\n', mat2str(size(raw_adc_data)));
fprintf('Frequencies shape: %s\n', mat2str(size(frequencies_hz)));

%% ===== EXACT PYTHON DATA CONVERSION =====
% Get user inputs (matching Python)
voltage_range_mv = 500;  % Default value
resolution_bits = 16;    % Default value
max_adc_value = (2^resolution_bits) / 2 - 1;

fprintf('Python-style conversion:\n');
fprintf('  Voltage range: %d mV\n', voltage_range_mv);
fprintf('  Resolution: %d bits\n', resolution_bits);
fprintf('  Max ADC value: %d\n', max_adc_value);

% EXACT Python conversion: data_mv = (raw_adc_data / max_adc_value) * voltage_range_mv
data_mv = (double(raw_adc_data) / max_adc_value) * voltage_range_mv;

fprintf('Converted data shape: %s\n', mat2str(size(data_mv)));
fprintf('Data range: [%.2f, %.2f] mV\n', min(data_mv(:)), max(data_mv(:)));

%% ===== EXACT PYTHON TIMING =====
% Python timing parameters
sample_interval_s = 24e-9;  % 24 ns
fs = 1 / sample_interval_s;
total_samples = size(data_mv, 2);
pre_trigger_samples = round(0.1 * total_samples);
time_ms = ((1:total_samples) - pre_trigger_samples) * (sample_interval_s * 1000);

fprintf('Timing parameters:\n');
fprintf('  Sample interval: %.2e s\n', sample_interval_s);
fprintf('  Sampling frequency: %.2e Hz\n', fs);
fprintf('  Total samples: %d\n', total_samples);
fprintf('  Pre-trigger samples: %d\n', pre_trigger_samples);

%% ===== EXACT PYTHON WINDOWS =====
% Python SNR windows
ECHO_WINDOW_MS = [0.35, 0.55];
NOISE_WINDOW_MS = [4.0, 5.0];

% Convert to sample indices (exactly like Python)
echo_start_sample = round((ECHO_WINDOW_MS(1) / 1000) * fs) + pre_trigger_samples;
echo_end_sample = round((ECHO_WINDOW_MS(2) / 1000) * fs) + pre_trigger_samples;
noise_start_sample = round((NOISE_WINDOW_MS(1) / 1000) * fs) + pre_trigger_samples;
noise_end_sample = round((NOISE_WINDOW_MS(2) / 1000) * fs) + pre_trigger_samples;

fprintf('Window calculations:\n');
fprintf('  Echo window: %.1f-%.1f ms (samples %d-%d)\n', ECHO_WINDOW_MS(1), ECHO_WINDOW_MS(2), echo_start_sample, echo_end_sample);
fprintf('  Noise window: %.1f-%.1f ms (samples %d-%d)\n', NOISE_WINDOW_MS(1), NOISE_WINDOW_MS(2), noise_start_sample, noise_end_sample);

%% ===== EXACT PYTHON FILTERING =====
% Python filter parameters
FILTER_CUTOFF_HZ = 60000;
FILTER_ORDER = 4;

% Apply exact Python filter
nyquist = 0.5 * fs;
normal_cutoff = FILTER_CUTOFF_HZ / nyquist;
[b, a] = butter(FILTER_ORDER, normal_cutoff, 'low');

fprintf('Filter parameters:\n');
fprintf('  Cutoff: %d Hz\n', FILTER_CUTOFF_HZ);
fprintf('  Order: %d\n', FILTER_ORDER);
fprintf('  Normalized cutoff: %.4f\n', normal_cutoff);

%% ===== EXACT PYTHON SNR CALCULATION =====
% Process all acquisitions (like Python)
all_snr_values = zeros(size(data_mv, 1), 1);

fprintf('\nCalculating SNR for all acquisitions...\n');

for i = 1:size(data_mv, 1)
    % Apply filter (exactly like Python)
    filtered_mv = filtfilt(b, a, data_mv(i, :));
    
    % Extract windows
    signal_segment = filtered_mv(echo_start_sample:echo_end_sample);
    noise_segment = filtered_mv(noise_start_sample:noise_end_sample);
    
    % Calculate RMS (exactly like Python)
    rms_signal = sqrt(mean(signal_segment.^2));
    rms_noise = sqrt(mean(noise_segment.^2));
    
    % Calculate SNR (exactly like Python)
    if rms_noise > 0
        snr_db = 20 * log10(rms_signal / rms_noise);
    else
        snr_db = Inf;
    end
    
    all_snr_values(i) = snr_db;
    
    if mod(i, 10) == 0 || i == 1
        fprintf('  Acquisition %d: SNR = %.2f dB\n', i, snr_db);
    end
end

%% ===== COMPARISON WITH PREVIOUS MATLAB ANALYSIS =====
fprintf('\n--- COMPARISON ANALYSIS ---\n');

% Load our previous analysis data
previous_data = load('/Users/deepshikhakaul/Documents/MATLAB/Optimized_071925_005/optimized_reconstruction_results.mat');

fprintf('Python-style analysis:\n');
fprintf('  Mean SNR: %.2f dB\n', mean(all_snr_values));
fprintf('  Min SNR: %.2f dB\n', min(all_snr_values));
fprintf('  Max SNR: %.2f dB\n', max(all_snr_values));
fprintf('  Std SNR: %.2f dB\n', std(all_snr_values));

% Check if we have SNR analysis from previous run
if isfield(previous_data, 'snr_analysis')
    fprintf('\nPrevious MATLAB analysis:\n');
    fprintf('  SNR: %.2f dB\n', previous_data.snr_analysis.original_snr_db);
else
    fprintf('\nPrevious MATLAB analysis: -13.2 dB (from console output)\n');
end

%% ===== DETAILED ANALYSIS =====
fprintf('\n--- DETAILED ANALYSIS ---\n');

% Analyze first acquisition in detail
acq_idx = 1;
filtered_mv = filtfilt(b, a, data_mv(acq_idx, :));
signal_segment = filtered_mv(echo_start_sample:echo_end_sample);
noise_segment = filtered_mv(noise_start_sample:noise_end_sample);

fprintf('Detailed analysis of acquisition %d:\n', acq_idx);
fprintf('  Signal segment length: %d samples\n', length(signal_segment));
fprintf('  Noise segment length: %d samples\n', length(noise_segment));
fprintf('  Signal RMS: %.6f mV\n', sqrt(mean(signal_segment.^2)));
fprintf('  Noise RMS: %.6f mV\n', sqrt(mean(noise_segment.^2)));
fprintf('  SNR: %.2f dB\n', all_snr_values(acq_idx));

% Plot comparison
figure(1);
subplot(2, 1, 1);
plot(time_ms, data_mv(acq_idx, :), 'b-', 'LineWidth', 1);
hold on;
plot(time_ms, filtered_mv, 'r-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
title('Raw vs Filtered Data (Acquisition 1)');
legend('Raw', 'Filtered');
grid on;

% Highlight windows
xlim([0, 6]);
ylim_vals = ylim;
patch([ECHO_WINDOW_MS(1), ECHO_WINDOW_MS(2), ECHO_WINDOW_MS(2), ECHO_WINDOW_MS(1)], ...
      [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
patch([NOISE_WINDOW_MS(1), NOISE_WINDOW_MS(2), NOISE_WINDOW_MS(2), NOISE_WINDOW_MS(1)], ...
      [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

subplot(2, 1, 2);
plot(1:length(all_snr_values), all_snr_values, 'bo-', 'LineWidth', 2);
xlabel('Acquisition Number');
ylabel('SNR (dB)');
title('SNR vs Acquisition Number (Python-style)');
grid on;

sgtitle('Python-style Analysis Results');
set(gcf, 'Color', 'w');

% Save results
save('PythonStyleAnalysis.mat', 'all_snr_values', 'data_mv', 'filtered_mv', 'signal_segment', 'noise_segment');

fprintf('\n--- ANALYSIS COMPLETE ---\n');
fprintf('Results saved to: PythonStyleAnalysis.mat\n'); 