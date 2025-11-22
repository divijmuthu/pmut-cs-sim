% Debug Reconstruction Pipeline
% This script checks each step of the reconstruction pipeline to find the issue
clearvars;
clc;
close all;

%% Load the latest results and check each component
fprintf('--- Debugging Reconstruction Pipeline ---\n');

% Load H matrix
h_matrix_file = '/Users/deepshikhakaul/Documents/MATLAB/Optimized_071925_010/optimized_H_matrix.mat';
if exist(h_matrix_file, 'file')
    fprintf('Loading H matrix from: %s\n', h_matrix_file);
    h_data = load(h_matrix_file);
    
    if isfield(h_data, 'H')
        H = h_data.H;
        fprintf('H matrix loaded: %s\n', mat2str(size(H)));
        fprintf('H matrix statistics:\n');
        fprintf('  Mean: %.6f\n', mean(H(:)));
        fprintf('  Std: %.6f\n', std(H(:)));
        fprintf('  Min: %.6f\n', min(H(:)));
        fprintf('  Max: %.6f\n', max(H(:)));
        fprintf('  Non-zero elements: %d/%d (%.2f%%)\n', sum(H(:) ~= 0), numel(H), 100*sum(H(:) ~= 0)/numel(H));
    else
        error('No H matrix found in file');
    end
else
    error('H matrix file not found');
end

% Load reconstruction results
reconstruction_file = '/Users/deepshikhakaul/Documents/MATLAB/Optimized_071925_010/optimized_reconstruction_results.mat';
if exist(reconstruction_file, 'file')
    fprintf('\nLoading reconstruction results from: %s\n', reconstruction_file);
    recon_data = load(reconstruction_file);
    
    if isfield(recon_data, 'image_recon')
        image_recon = recon_data.image_recon;
        fprintf('Reconstruction image: %s\n', mat2str(size(image_recon)));
        fprintf('Image statistics:\n');
        fprintf('  Mean: %.6f\n', mean(image_recon(:)));
        fprintf('  Std: %.6f\n', std(image_recon(:)));
        fprintf('  Min: %.6f\n', min(image_recon(:)));
        fprintf('  Max: %.6f\n', max(image_recon(:)));
        fprintf('  Non-zero pixels: %d/%d (%.2f%%)\n', sum(image_recon(:) ~= 0), numel(image_recon), 100*sum(image_recon(:) ~= 0)/numel(image_recon));
    end
    
    if isfield(recon_data, 'convergence_history')
        convergence_history = recon_data.convergence_history;
        fprintf('Convergence history: %d iterations\n', length(convergence_history));
        if ~isempty(convergence_history)
            fprintf('  Final residual: %.2e\n', convergence_history(end));
        end
    end
else
    error('Reconstruction file not found');
end

%% Check if the issue is in the measurement vector preparation
fprintf('\n--- Checking Measurement Vector Preparation ---\n');

% Load the original data to recreate the measurement vector
h5_file_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/data_2025-07-18_15-08-13_subtracted.h5';
raw_adc_data = h5read(h5_file_path, '/echo_data_raw_adc');
if size(raw_adc_data, 1) > size(raw_adc_data, 2)
    raw_adc_data = raw_adc_data';
end

% Convert to voltage
voltage_range_mv = 500;
resolution_bits = 16;
max_adc_value = (2^resolution_bits) / 2 - 1;
data_mv = (double(raw_adc_data) / max_adc_value) * voltage_range_mv;

% Use last 100 acquisitions
last_100_data = data_mv(101:200, :);

% Apply preprocessing (same as main script)
echo_data = last_100_data;
for acq = 1:size(echo_data, 1)
    signal = echo_data(acq, :);
    baseline = mean(signal(1:min(1000, length(signal))));
    echo_data(acq, :) = signal - baseline;
end

% Apply lowpass filter
fs = 41.67e6;
cutoff_freq = 60000;
nyquist_freq = fs / 2;
normalized_cutoff = cutoff_freq / nyquist_freq;
[b, a] = butter(4, normalized_cutoff, 'low');

echo_data_filtered = zeros(size(echo_data));
for acq = 1:size(echo_data, 1)
    echo_data_filtered(acq, :) = filtfilt(b, a, echo_data(acq, :));
end

% Apply SNR scaling (1.67x from the output)
scaling_factor = 1.67;
echo_data_scaled = echo_data_filtered * scaling_factor;

% Downsample
downsample_factor = 10;
echo_data_downsampled = zeros(size(echo_data_scaled, 1), floor(size(echo_data_scaled, 2) / downsample_factor));
for acq = 1:size(echo_data_scaled, 1)
    echo_data_downsampled(acq, :) = decimate(echo_data_scaled(acq, :), downsample_factor);
end

% Normalize
for acq = 1:size(echo_data_downsampled, 1)
    signal = echo_data_downsampled(acq, :);
    signal_std = std(signal);
    if signal_std > 0
        echo_data_downsampled(acq, :) = signal / signal_std;
    end
end

fprintf('Preprocessed data statistics:\n');
fprintf('  Mean: %.6f\n', mean(echo_data_downsampled(:)));
fprintf('  Std: %.6f\n', std(echo_data_downsampled(:)));
fprintf('  Min: %.6f\n', min(echo_data_downsampled(:)));
fprintf('  Max: %.6f\n', max(echo_data_downsampled(:)));

%% Test feature extraction
fprintf('\n--- Testing Feature Extraction ---\n');

% Extract features (same as main script)
features = zeros(size(echo_data_downsampled, 1), 1);
for acq = 1:size(echo_data_downsampled, 1)
    signal = echo_data_downsampled(acq, :);
    
    % Method 1: Energy-based feature
    signal_energy = sum(signal.^2);
    
    % Method 2: Peak detection
    [peaks, peak_locs] = findpeaks(abs(signal), 'MinPeakHeight', 0.1*max(abs(signal)));
    if ~isempty(peaks)
        peak_energy = sum(peaks.^2);
    else
        peak_energy = 0;
    end
    
    % Method 3: Cross-correlation
    t = 1:length(signal);
    expected_pattern = exp(-((t - length(signal)/2).^2) / (2 * (length(signal)/8)^2));
    correlation = xcorr(signal, expected_pattern);
    correlation_peak = max(abs(correlation));
    
    % Method 4: Frequency domain
    signal_fft = abs(fft(signal));
    dominant_freq_energy = max(signal_fft(1:end/2));
    
    % Combine features
    feature_value = 0.3 * signal_energy + 0.3 * peak_energy + 0.2 * correlation_peak + 0.2 * dominant_freq_energy;
    features(acq) = abs(feature_value);
end

fprintf('Feature extraction results:\n');
fprintf('  Mean feature: %.6f\n', mean(features));
fprintf('  Std feature: %.6f\n', std(features));
fprintf('  Min feature: %.6f\n', min(features));
fprintf('  Max feature: %.6f\n', max(features));
fprintf('  Non-zero features: %d/%d (%.2f%%)\n', sum(features ~= 0), length(features), 100*sum(features ~= 0)/length(features));

% Normalize features
if norm(features) > 0
    features_normalized = features / norm(features);
else
    features_normalized = features;
end

fprintf('Normalized features:\n');
fprintf('  Mean: %.6f\n', mean(features_normalized));
fprintf('  Std: %.6f\n', std(features_normalized));
fprintf('  Min: %.6f\n', min(features_normalized));
fprintf('  Max: %.6f\n', max(features_normalized));

%% Test ADMM reconstruction with the extracted features
fprintf('\n--- Testing ADMM Reconstruction ---\n');

% Use the extracted features as measurement vector
y = features_normalized;

fprintf('Measurement vector statistics:\n');
fprintf('  Size: %s\n', mat2str(size(y)));
fprintf('  Mean: %.6f\n', mean(y));
fprintf('  Std: %.6f\n', std(y));
fprintf('  Min: %.6f\n', min(y));
fprintf('  Max: %.6f\n', max(y));
fprintf('  Norm: %.6f\n', norm(y));

% Check H matrix condition
cond_num = cond(H);
fprintf('H matrix condition number: %.2e\n', cond_num);

% Check if H*y is reasonable
Hy = H * y;
fprintf('H*y statistics:\n');
fprintf('  Mean: %.6f\n', mean(Hy));
fprintf('  Std: %.6f\n', std(Hy));
fprintf('  Min: %.6f\n', min(Hy));
fprintf('  Max: %.6f\n', max(Hy));

% Test simple least squares
fprintf('\n--- Testing Simple Least Squares ---\n');
x_ls = H \ y;
fprintf('Least squares solution:\n');
fprintf('  Mean: %.6f\n', mean(x_ls));
fprintf('  Std: %.6f\n', std(x_ls));
fprintf('  Min: %.6f\n', min(x_ls));
fprintf('  Max: %.6f\n', max(x_ls));
fprintf('  Non-zero elements: %d/%d (%.2f%%)\n', sum(x_ls ~= 0), length(x_ls), 100*sum(x_ls ~= 0)/length(x_ls));

% Reshape to image
image_ls = reshape(x_ls, 51, 76);
fprintf('Least squares image:\n');
fprintf('  Mean: %.6f\n', mean(image_ls(:)));
fprintf('  Std: %.6f\n', std(image_ls(:)));
fprintf('  Min: %.6f\n', min(image_ls(:)));
fprintf('  Max: %.6f\n', max(image_ls(:)));

%% Plot results
figure(1);
set(gcf, 'Position', [100, 100, 1400, 800]);

% Original reconstruction
subplot(2, 4, 1);
imagesc(image_recon);
colorbar;
title('Original Reconstruction (All Zeros)');
axis equal;

% Least squares reconstruction
subplot(2, 4, 2);
imagesc(image_ls);
colorbar;
title('Least Squares Reconstruction');
axis equal;

% H matrix
subplot(2, 4, 3);
imagesc(H);
colorbar;
title('H Matrix');
xlabel('Grid Points'); ylabel('Acquisitions');

% Measurement vector
subplot(2, 4, 4);
plot(y, 'b-', 'LineWidth', 2);
title('Measurement Vector');
xlabel('Acquisition'); ylabel('Feature Value');
grid on;

% Feature extraction
subplot(2, 4, 5);
plot(features, 'r-', 'LineWidth', 2);
title('Extracted Features');
xlabel('Acquisition'); ylabel('Feature Value');
grid on;

% H*y
subplot(2, 4, 6);
plot(Hy, 'g-', 'LineWidth', 2);
title('H*y');
xlabel('Grid Point'); ylabel('Value');
grid on;

% Histogram of H matrix
subplot(2, 4, 7);
histogram(H(:), 50);
title('H Matrix Distribution');
xlabel('Value'); ylabel('Count');
grid on;

% Histogram of measurement vector
subplot(2, 4, 8);
histogram(y, 20);
title('Measurement Vector Distribution');
xlabel('Value'); ylabel('Count');
grid on;

sgtitle('Reconstruction Pipeline Debug');
set(gcf, 'Color', 'w');

fprintf('\n--- Debug Complete ---\n');
fprintf('Check the plots to see where the pipeline is failing!\n'); 