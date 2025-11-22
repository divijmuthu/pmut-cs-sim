% Debug SNR Calculation
% This script debugs why the SNR calculation is still wrong
clearvars;
clc;
close all;

%% Load the same data as our working analysis
fprintf('--- Debugging SNR Calculation ---\n');

% Load data (same as CheckLast100SNR.m)
h5_file_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/data_2025-07-18_15-08-13_subtracted.h5';
raw_adc_data = h5read(h5_file_path, '/echo_data_raw_adc');
if size(raw_adc_data, 1) > size(raw_adc_data, 2)
    raw_adc_data = raw_adc_data';
end

% Convert to voltage (same as working analysis)
voltage_range_mv = 500;
resolution_bits = 16;
max_adc_value = (2^resolution_bits) / 2 - 1;
data_mv = (double(raw_adc_data) / max_adc_value) * voltage_range_mv;

% Use last 100 acquisitions
last_100_data = data_mv(101:200, :);
fprintf('Analyzing last 100 acquisitions (indices 101-200)\n');

%% Apply the same preprocessing as the main script
fprintf('\n--- Applying Main Script Preprocessing ---\n');

% Step 1: Baseline correction
echo_data = last_100_data;
for acq = 1:size(echo_data, 1)
    signal = echo_data(acq, :);
    baseline = mean(signal(1:min(1000, length(signal))));
    echo_data(acq, :) = signal - baseline;
end
fprintf('Step 1: Baseline correction completed.\n');

% Step 2: Lowpass filter
fs = 41.67e6;  % Original sampling rate
cutoff_freq = 60000;
nyquist_freq = fs / 2;
normalized_cutoff = cutoff_freq / nyquist_freq;
[b, a] = butter(4, normalized_cutoff, 'low');

echo_data_filtered = zeros(size(echo_data));
for acq = 1:size(echo_data, 1)
    echo_data_filtered(acq, :) = filtfilt(b, a, echo_data(acq, :));
end
fprintf('Step 2: Lowpass filtering completed.\n');

%% Now calculate SNR using the main script's method
fprintf('\n--- Calculating SNR (Main Script Method) ---\n');

% Timing parameters (same as main script)
sample_interval_s = 24e-9;
fs_main = 1 / sample_interval_s;
total_samples = size(echo_data_filtered, 2);
pre_trigger_samples = round(0.1 * total_samples);

% Windows (same as main script)
ECHO_WINDOW_MS = [0.35, 0.55];
NOISE_WINDOW_MS = [4.0, 5.0];

echo_start_sample = round((ECHO_WINDOW_MS(1) / 1000) * fs_main) + pre_trigger_samples;
echo_end_sample = round((ECHO_WINDOW_MS(2) / 1000) * fs_main) + pre_trigger_samples;
noise_start_sample = round((NOISE_WINDOW_MS(1) / 1000) * fs_main) + pre_trigger_samples;
noise_end_sample = round((NOISE_WINDOW_MS(2) / 1000) * fs_main) + pre_trigger_samples;

fprintf('Main script timing:\n');
fprintf('  Echo window: %.1f-%.1f ms (samples %d-%d)\n', ECHO_WINDOW_MS(1), ECHO_WINDOW_MS(2), echo_start_sample, echo_end_sample);
fprintf('  Noise window: %.1f-%.1f ms (samples %d-%d)\n', NOISE_WINDOW_MS(1), NOISE_WINDOW_MS(2), noise_start_sample, noise_end_sample);

% Calculate SNR for first acquisition
acq_idx = 1;
signal_segment = echo_data_filtered(acq_idx, echo_start_sample:echo_end_sample);
noise_segment = echo_data_filtered(acq_idx, noise_start_sample:noise_end_sample);

echo_power = mean(signal_segment.^2);
noise_power = mean(noise_segment.^2);

if noise_power > 0
    snr_db = 10 * log10(echo_power / noise_power);
else
    snr_db = Inf;
end

fprintf('Main script SNR calculation:\n');
fprintf('  Echo power: %.2e\n', echo_power);
fprintf('  Noise power: %.2e\n', noise_power);
fprintf('  SNR: %.2f dB\n', snr_db);

%% Compare with our working analysis
fprintf('\n--- Comparing with Working Analysis ---\n');

% Working analysis timing (from CheckLast100SNR.m)
fs_working = 41.67e6;
sample_interval_s_working = 1 / fs_working;
pre_trigger_samples_working = round(0.1 * size(last_100_data, 2));

echo_start_sample_working = round((ECHO_WINDOW_MS(1) / 1000) * fs_working) + pre_trigger_samples_working;
echo_end_sample_working = round((ECHO_WINDOW_MS(2) / 1000) * fs_working) + pre_trigger_samples_working;
noise_start_sample_working = round((NOISE_WINDOW_MS(1) / 1000) * fs_working) + pre_trigger_samples_working;
noise_end_sample_working = round((NOISE_WINDOW_MS(2) / 1000) * fs_working) + pre_trigger_samples_working;

fprintf('Working analysis timing:\n');
fprintf('  Echo window: %.1f-%.1f ms (samples %d-%d)\n', ECHO_WINDOW_MS(1), ECHO_WINDOW_MS(2), echo_start_sample_working, echo_end_sample_working);
fprintf('  Noise window: %.1f-%.1f ms (samples %d-%d)\n', NOISE_WINDOW_MS(1), NOISE_WINDOW_MS(2), noise_start_sample_working, noise_end_sample_working);

% Apply filter to raw data (same as working analysis)
echo_data_filtered_working = zeros(size(last_100_data));
for acq = 1:size(last_100_data, 1)
    echo_data_filtered_working(acq, :) = filtfilt(b, a, last_100_data(acq, :));
end

% Calculate SNR using working method
signal_segment_working = echo_data_filtered_working(acq_idx, echo_start_sample_working:echo_end_sample_working);
noise_segment_working = echo_data_filtered_working(acq_idx, noise_start_sample_working:noise_end_sample_working);

echo_power_working = mean(signal_segment_working.^2);
noise_power_working = mean(noise_segment_working.^2);

if noise_power_working > 0
    snr_db_working = 10 * log10(echo_power_working / noise_power_working);
else
    snr_db_working = Inf;
end

fprintf('Working analysis SNR calculation:\n');
fprintf('  Echo power: %.2e\n', echo_power_working);
fprintf('  Noise power: %.2e\n', noise_power_working);
fprintf('  SNR: %.2f dB\n', snr_db_working);

%% Plot comparison
figure(1);
subplot(2, 1, 1);
plot(echo_data_filtered(acq_idx, :), 'b-', 'LineWidth', 1);
hold on;
plot(echo_data_filtered_working(acq_idx, :), 'r-', 'LineWidth', 1);
xlabel('Sample Index');
ylabel('Amplitude (mV)');
title('Filtered Data Comparison (Acquisition 1)');
legend('Main Script', 'Working Analysis');
grid on;

% Highlight windows
ylim_vals = ylim;
patch([echo_start_sample, echo_end_sample, echo_end_sample, echo_start_sample], ...
      [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
patch([noise_start_sample, noise_end_sample, noise_end_sample, noise_start_sample], ...
      [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

subplot(2, 1, 2);
plot(echo_start_sample:echo_end_sample, signal_segment, 'b-', 'LineWidth', 2);
hold on;
plot(noise_start_sample:noise_end_sample, noise_segment, 'r-', 'LineWidth', 2);
xlabel('Sample Index');
ylabel('Amplitude (mV)');
title('Signal vs Noise Segments (Main Script)');
legend('Signal', 'Noise');
grid on;

sgtitle('SNR Calculation Debug');
set(gcf, 'Color', 'w');

fprintf('\n--- Debug Complete ---\n');
fprintf('The issue is likely in the timing window calculations!\n'); 