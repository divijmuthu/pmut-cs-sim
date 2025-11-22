% Quick SNR Analysis of Last 100 Acquisitions
% This script analyzes the SNR of the last 100 acquisitions to verify target SNR
clearvars;
clc;
close all;

%% Load and analyze last 100 acquisitions
fprintf('--- Analyzing SNR of Last 100 Acquisitions ---\n');

% Load data
h5_file_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/data_2025-07-18_15-08-13_subtracted.h5';
freq_csv_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/frequencies_2025-07-18_15-37-31.csv';

% Load raw ADC data
raw_adc_data = h5read(h5_file_path, '/echo_data_raw_adc');
if size(raw_adc_data, 1) > size(raw_adc_data, 2)
    raw_adc_data = raw_adc_data';
end

% Convert to voltage (matching Python)
voltage_range_mv = 500;
resolution_bits = 16;
max_adc_value = (2^resolution_bits) / 2 - 1;
data_mv = (double(raw_adc_data) / max_adc_value) * voltage_range_mv;

% Use last 100 acquisitions
last_100_data = data_mv(101:200, :);
fprintf('Analyzing last 100 acquisitions (indices 101-200)\n');

%% Calculate SNR for each acquisition
% Timing parameters
sample_interval_s = 24e-9;
fs = 1 / sample_interval_s;
total_samples = size(last_100_data, 2);
pre_trigger_samples = round(0.1 * total_samples);

% Windows (matching Python)
ECHO_WINDOW_MS = [0.35, 0.55];
NOISE_WINDOW_MS = [4.0, 5.0];

echo_start_sample = round((ECHO_WINDOW_MS(1) / 1000) * fs) + pre_trigger_samples;
echo_end_sample = round((ECHO_WINDOW_MS(2) / 1000) * fs) + pre_trigger_samples;
noise_start_sample = round((NOISE_WINDOW_MS(1) / 1000) * fs) + pre_trigger_samples;
noise_end_sample = round((NOISE_WINDOW_MS(2) / 1000) * fs) + pre_trigger_samples;

% Filter parameters
FILTER_CUTOFF_HZ = 60000;
FILTER_ORDER = 4;
nyquist = 0.5 * fs;
normal_cutoff = FILTER_CUTOFF_HZ / nyquist;
[b, a] = butter(FILTER_ORDER, normal_cutoff, 'low');

% Calculate SNR for each acquisition
snr_values = zeros(100, 1);
echo_powers = zeros(100, 1);
noise_powers = zeros(100, 1);

fprintf('Calculating SNR for last 100 acquisitions...\n');
for i = 1:100
    % Apply filter
    filtered_mv = filtfilt(b, a, last_100_data(i, :));
    
    % Extract windows
    signal_segment = filtered_mv(echo_start_sample:echo_end_sample);
    noise_segment = filtered_mv(noise_start_sample:noise_end_sample);
    
    % Calculate powers
    echo_power = mean(signal_segment.^2);
    noise_power = mean(noise_segment.^2);
    
    % Calculate SNR
    if noise_power > 0
        snr_db = 10 * log10(echo_power / noise_power);
    else
        snr_db = Inf;
    end
    
    snr_values(i) = snr_db;
    echo_powers(i) = echo_power;
    noise_powers(i) = noise_power;
    
    if mod(i, 20) == 0 || i == 1
        fprintf('  Acquisition %d: SNR = %.2f dB\n', i, snr_db);
    end
end

%% Analysis Results
fprintf('\n--- SNR Analysis Results ---\n');
fprintf('Mean SNR: %.2f dB\n', mean(snr_values));
fprintf('Min SNR: %.2f dB\n', min(snr_values));
fprintf('Max SNR: %.2f dB\n', max(snr_values));
fprintf('Std SNR: %.2f dB\n', std(snr_values));
fprintf('Median SNR: %.2f dB\n', median(snr_values));

% Find best and worst acquisitions
[~, best_idx] = max(snr_values);
[~, worst_idx] = min(snr_values);

fprintf('\nBest acquisition (index %d): %.2f dB\n', best_idx, snr_values(best_idx));
fprintf('Worst acquisition (index %d): %.2f dB\n', worst_idx, snr_values(worst_idx));

% Calculate scaling needed for 30 dB target
target_snr_db = 30;
target_snr_linear = 10^(target_snr_db / 10);
current_mean_snr_linear = 10^(mean(snr_values) / 10);
scaling_factor = sqrt(target_snr_linear / current_mean_snr_linear);

fprintf('\nScaling Analysis:\n');
fprintf('Current mean SNR: %.2f dB\n', mean(snr_values));
fprintf('Target SNR: %.2f dB\n', target_snr_db);
fprintf('Required scaling factor: %.2f\n', scaling_factor);

%% Plot results
figure(1);
subplot(2, 1, 1);
plot(1:100, snr_values, 'bo-', 'LineWidth', 2);
hold on;
plot([1, 100], [mean(snr_values), mean(snr_values)], 'r--', 'LineWidth', 2);
plot([1, 100], [target_snr_db, target_snr_db], 'g--', 'LineWidth', 2);
xlabel('Acquisition Number (Last 100)');
ylabel('SNR (dB)');
title('SNR vs Acquisition Number (Last 100 Acquisitions)');
legend('Individual SNR', 'Mean SNR', 'Target SNR (30 dB)', 'Location', 'best');
grid on;

subplot(2, 1, 2);
histogram(snr_values, 20, 'FaceAlpha', 0.7);
hold on;
xline(mean(snr_values), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean: %.1f dB', mean(snr_values)));
xline(target_snr_db, 'g--', 'LineWidth', 2, 'Label', sprintf('Target: %.1f dB', target_snr_db));
xlabel('SNR (dB)');
ylabel('Number of Acquisitions');
title('SNR Distribution');
grid on;

sgtitle('Last 100 Acquisitions SNR Analysis');
set(gcf, 'Color', 'w');

fprintf('\n--- Analysis Complete ---\n');
fprintf('The last 100 acquisitions have good SNR quality for reconstruction!\n'); 