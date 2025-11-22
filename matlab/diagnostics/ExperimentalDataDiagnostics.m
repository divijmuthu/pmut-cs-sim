% Experimental Data Diagnostics - Analyze SNR and Data Quality Issues
clearvars;
clc;
close all;

%% ===== LOAD AND ANALYZE EXPERIMENTAL DATA =====
fprintf('--- Experimental Data Diagnostics ---\n');

% Load experimental data
h5_file_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/data_2025-07-18_15-08-13_subtracted.h5';
delay_csv_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/delays_2025-07-18_15-37-31.csv';
freq_csv_path = '/Users/deepshikhakaul/Documents/pmut-cs-sim/frequencies_2025-07-18_15-37-31.csv';

% Load raw data
echo_data_raw = h5read(h5_file_path, '/echo_data_raw_adc');
if size(echo_data_raw, 1) > size(echo_data_raw, 2)
    echo_data_raw = echo_data_raw';
end

% Load delay and frequency data
delay_profiles = table2array(readtable(delay_csv_path));
frequencies = table2array(readtable(freq_csv_path));

% Use last 100 acquisitions
echo_data = echo_data_raw(end-99:end, :);
delay_profiles = delay_profiles(end-99:end, :);
frequencies = frequencies(end-99:end, :);

fprintf('Data loaded: %d acquisitions x %d samples\n', size(echo_data, 1), size(echo_data, 2));

%% ===== COMPREHENSIVE DATA ANALYSIS =====

% 1. Raw Data Analysis
fprintf('\n--- Raw Data Analysis ---\n');
analyze_raw_data(echo_data);

% 2. Preprocessing Analysis
fprintf('\n--- Preprocessing Analysis ---\n');
[echo_data_processed, preprocessing_stats] = analyze_preprocessing(echo_data);

% 3. SNR Analysis
fprintf('\n--- SNR Analysis ---\n');
snr_analysis = analyze_snr_characteristics(echo_data_processed);

% 4. Comparison with Simulation Standards
fprintf('\n--- Comparison with Simulation Standards ---\n');
compare_with_simulation_standards(snr_analysis);

% 5. Generate Diagnostic Plots
fprintf('\n--- Generating Diagnostic Plots ---\n');
generate_diagnostic_plots(echo_data, echo_data_processed, snr_analysis, preprocessing_stats);

fprintf('\n--- Diagnostics Complete ---\n');

%% ===== HELPER FUNCTIONS =====

function analyze_raw_data(echo_data)
    % Analyze raw experimental data characteristics
    
    fprintf('Raw Data Statistics:\n');
    fprintf('  Data Type: %s\n', class(echo_data));
    fprintf('  Data Range: [%.2e, %.2e]\n', min(echo_data(:)), max(echo_data(:)));
    fprintf('  Mean Value: %.2e\n', mean(echo_data(:)));
    fprintf('  Standard Deviation: %.2e\n', std(echo_data(:)));
    fprintf('  Dynamic Range: %.1f dB\n', 20*log10(max(abs(echo_data(:))) / (min(abs(echo_data(:))) + eps)));
    
    % Check for data quality issues
    zero_count = sum(echo_data(:) == 0);
    total_elements = numel(echo_data);
    zero_percentage = 100 * zero_count / total_elements;
    
    fprintf('  Zero Values: %d/%d (%.1f%%)\n', zero_count, total_elements, zero_percentage);
    
    % Check for saturation
    unique_values = unique(echo_data(:));
    fprintf('  Unique Values: %d\n', length(unique_values));
    
    if length(unique_values) < 100
        warning('Low number of unique values suggests possible saturation or quantization issues');
    end
end

function [echo_data_processed, stats] = analyze_preprocessing(echo_data)
    % Analyze preprocessing steps and their effects
    
    stats = struct();
    
    % Step 1: Convert to double
    echo_data_double = double(echo_data);
    stats.original_range = [min(echo_data(:)), max(echo_data(:))];
    stats.double_range = [min(echo_data_double(:)), max(echo_data_double(:))];
    
    % Step 2: Baseline correction
    echo_data_baseline = echo_data_double;
    for acq = 1:size(echo_data_baseline, 1)
        signal = echo_data_baseline(acq, :);
        baseline = mean(signal(1:min(1000, length(signal))));
        echo_data_baseline(acq, :) = signal - baseline;
    end
    stats.baseline_corrected_range = [min(echo_data_baseline(:)), max(echo_data_baseline(:))];
    
    % Step 3: Lowpass filter (60 kHz)
    fs = 2e6;  % 2 MHz sampling
    cutoff = 60e3;  % 60 kHz
    [b, a] = butter(4, cutoff/(fs/2), 'low');
    
    echo_data_filtered = zeros(size(echo_data_baseline));
    for acq = 1:size(echo_data_baseline, 1)
        echo_data_filtered(acq, :) = filtfilt(b, a, echo_data_baseline(acq, :));
    end
    stats.filtered_range = [min(echo_data_filtered(:)), max(echo_data_filtered(:))];
    
    % Step 4: Downsampling
    downsample_factor = 10;
    echo_data_downsampled = zeros(size(echo_data_filtered, 1), floor(size(echo_data_filtered, 2)/downsample_factor));
    for acq = 1:size(echo_data_filtered, 1)
        echo_data_downsampled(acq, :) = decimate(echo_data_filtered(acq, :), downsample_factor);
    end
    stats.downsampled_range = [min(echo_data_downsampled(:)), max(echo_data_downsampled(:))];
    
    % Step 5: Normalization
    echo_data_normalized = zeros(size(echo_data_downsampled));
    for acq = 1:size(echo_data_downsampled, 1)
        signal = echo_data_downsampled(acq, :);
        signal_std = std(signal);
        if signal_std > 0
            echo_data_normalized(acq, :) = signal / signal_std;
        else
            echo_data_normalized(acq, :) = signal;
        end
    end
    stats.normalized_range = [min(echo_data_normalized(:)), max(echo_data_normalized(:))];
    
    echo_data_processed = echo_data_normalized;
    
    fprintf('Preprocessing Effects:\n');
    fprintf('  Original Range: [%.2e, %.2e]\n', stats.original_range(1), stats.original_range(2));
    fprintf('  Double Range: [%.2e, %.2e]\n', stats.double_range(1), stats.double_range(2));
    fprintf('  Baseline Corrected: [%.2e, %.2e]\n', stats.baseline_corrected_range(1), stats.baseline_corrected_range(2));
    fprintf('  Filtered Range: [%.2e, %.2e]\n', stats.filtered_range(1), stats.filtered_range(2));
    fprintf('  Downsampled Range: [%.2e, %.2e]\n', stats.downsampled_range(1), stats.downsampled_range(2));
    fprintf('  Normalized Range: [%.2e, %.2e]\n', stats.normalized_range(1), stats.normalized_range(2));
end

function snr_analysis = analyze_snr_characteristics(echo_data)
    % Comprehensive SNR analysis
    
    snr_analysis = struct();
    
    % Calculate signal statistics
    signal_power = mean(echo_data(:).^2);
    signal_rms = rms(echo_data(:));
    signal_peak = max(abs(echo_data(:)));
    signal_dynamic_range = 20 * log10(signal_peak / (min(abs(echo_data(:))) + eps));
    
    % Estimate noise floor (using first 1000 samples as noise reference)
    noise_samples = echo_data(:, 1:min(1000, size(echo_data, 2)));
    noise_power = mean(noise_samples(:).^2);
    noise_rms = rms(noise_samples(:));
    
    % Calculate SNR
    if noise_power > 0
        snr_linear = signal_power / noise_power;
        snr_db = 10 * log10(snr_linear);
    else
        snr_db = Inf;
        snr_linear = Inf;
    end
    
    % Store results
    snr_analysis.signal_power = signal_power;
    snr_analysis.signal_rms = signal_rms;
    snr_analysis.signal_peak = signal_peak;
    snr_analysis.dynamic_range_db = signal_dynamic_range;
    snr_analysis.noise_power = noise_power;
    snr_analysis.noise_rms = noise_rms;
    snr_analysis.snr_db = snr_db;
    snr_analysis.snr_linear = snr_linear;
    
    fprintf('SNR Analysis:\n');
    fprintf('  Signal Power: %.2e\n', signal_power);
    fprintf('  Signal RMS: %.2e\n', signal_rms);
    fprintf('  Signal Peak: %.2e\n', signal_peak);
    fprintf('  Dynamic Range: %.1f dB\n', signal_dynamic_range);
    fprintf('  Noise Power: %.2e\n', noise_power);
    fprintf('  Noise RMS: %.2e\n', noise_rms);
    fprintf('  SNR: %.1f dB (%.2e linear)\n', snr_db, snr_linear);
    
    % Analyze signal quality
    if snr_db < 10
        fprintf('  ⚠️  WARNING: Very low SNR (< 10 dB)\n');
    elseif snr_db < 20
        fprintf('  ⚠️  WARNING: Low SNR (< 20 dB)\n');
    elseif snr_db < 30
        fprintf('  ⚠️  WARNING: Moderate SNR (< 30 dB)\n');
    else
        fprintf('  ✅ Good SNR (≥ 30 dB)\n');
    end
end

function compare_with_simulation_standards(snr_analysis)
    % Compare experimental SNR with simulation standards
    
    fprintf('Comparison with Simulation Standards:\n');
    
    % Simulation standards (from the codebase)
    simulation_snr_targets = [25, 30, 35, 40, 90];  % dB values found in simulations
    experimental_snr = snr_analysis.snr_db;
    
    fprintf('  Experimental SNR: %.1f dB\n', experimental_snr);
    fprintf('  Simulation Targets: %s dB\n', mat2str(simulation_snr_targets));
    
    % Find closest simulation target
    [~, closest_idx] = min(abs(simulation_snr_targets - experimental_snr));
    closest_target = simulation_snr_targets(closest_idx);
    
    fprintf('  Closest Simulation Target: %.1f dB\n', closest_target);
    
    % Calculate improvement needed
    improvement_needed = closest_target - experimental_snr;
    
    if improvement_needed > 0
        fprintf('  📈 SNR Improvement Needed: +%.1f dB\n', improvement_needed);
        
        % Suggest scaling factor
        scaling_factor = 10^(improvement_needed / 20);
        fprintf('  🔧 Suggested Scaling Factor: %.1fx\n', scaling_factor);
    else
        fprintf('  ✅ SNR meets simulation standards\n');
    end
    
    % Analyze root causes
    fprintf('\n  Root Cause Analysis:\n');
    if experimental_snr < 10
        fprintf('    • Very poor signal quality - possible hardware issues\n');
        fprintf('    • Check transducer connections and amplifier settings\n');
    elseif experimental_snr < 20
        fprintf('    • Low signal quality - possible gain/amplification issues\n');
        fprintf('    • Consider increasing amplifier gain or signal averaging\n');
    elseif experimental_snr < 30
        fprintf('    • Moderate signal quality - may need signal processing improvements\n');
        fprintf('    • Consider better filtering or noise reduction techniques\n');
    end
end

function generate_diagnostic_plots(echo_data_raw, echo_data_processed, snr_analysis, preprocessing_stats)
    % Generate comprehensive diagnostic plots
    
    figure(1);
    set(gcf, 'Position', [100, 100, 1200, 800]);
    clf;
    
    % Plot 1: Raw vs Processed Data
    subplot(2, 3, 1);
    sample_acq = 1;
    t_raw = (1:size(echo_data_raw, 2)) / 2e6 * 1e6;  % Time in microseconds
    t_processed = (1:size(echo_data_processed, 2)) / (2e6/10) * 1e6;  % After downsampling
    
    plot(t_raw(1:min(1000, length(t_raw))), echo_data_raw(sample_acq, 1:min(1000, size(echo_data_raw, 2))));
    hold on;
    plot(t_processed(1:min(1000, length(t_processed))), echo_data_processed(sample_acq, 1:min(1000, size(echo_data_processed, 2))));
    title('Raw vs Processed Data (Acquisition 1)');
    xlabel('Time (μs)'); ylabel('Amplitude');
    legend('Raw', 'Processed');
    grid on;
    
    % Plot 2: SNR Analysis
    subplot(2, 3, 2);
    snr_metrics = [snr_analysis.signal_power, snr_analysis.noise_power, snr_analysis.snr_db];
    bar(snr_metrics);
    set(gca, 'XTickLabel', {'Signal Power', 'Noise Power', 'SNR (dB)'});
    title('Signal Quality Metrics');
    ylabel('Value');
    grid on;
    
    % Plot 3: Data Distribution
    subplot(2, 3, 3);
    histogram(echo_data_processed(:), 50, 'Normalization', 'probability');
    title('Processed Data Distribution');
    xlabel('Amplitude'); ylabel('Probability');
    grid on;
    
    % Plot 4: Preprocessing Effects
    subplot(2, 3, 4);
    ranges = [preprocessing_stats.original_range; preprocessing_stats.double_range; ...
              preprocessing_stats.baseline_corrected_range; preprocessing_stats.filtered_range; ...
              preprocessing_stats.downsampled_range; preprocessing_stats.normalized_range];
    bar(ranges);
    set(gca, 'XTickLabel', {'Original', 'Double', 'Baseline', 'Filtered', 'Downsampled', 'Normalized'});
    title('Preprocessing Effects on Data Range');
    ylabel('Range');
    legend('Min', 'Max');
    grid on;
    
    % Plot 5: Signal Quality Over Acquisitions
    subplot(2, 3, 5);
    snr_per_acq = zeros(size(echo_data_processed, 1), 1);
    for acq = 1:size(echo_data_processed, 1)
        signal = echo_data_processed(acq, :);
        signal_power = mean(signal.^2);
        noise_power = mean(signal(1:min(1000, length(signal))).^2);
        if noise_power > 0
            snr_per_acq(acq) = 10 * log10(signal_power / noise_power);
        else
            snr_per_acq(acq) = 0;
        end
    end
    plot(1:length(snr_per_acq), snr_per_acq, 'o-');
    title('SNR per Acquisition');
    xlabel('Acquisition Number'); ylabel('SNR (dB)');
    grid on;
    
    % Plot 6: Summary
    subplot(2, 3, 6);
    summary_data = [snr_analysis.signal_peak, snr_analysis.noise_rms, snr_analysis.dynamic_range_db];
    bar(summary_data);
    set(gca, 'XTickLabel', {'Peak Signal', 'Noise RMS', 'Dynamic Range (dB)'});
    title('Data Quality Summary');
    ylabel('Value');
    grid on;
    
    sgtitle('Experimental Data Diagnostics');
    set(gcf, 'Color', 'w');
    
    % Save diagnostic plot
    saveas(gcf, 'ExperimentalDataDiagnostics.png');
    fprintf('  Diagnostic plots saved to: ExperimentalDataDiagnostics.png\n');
end 