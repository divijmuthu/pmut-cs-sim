% =========================================================================
% EXPERIMENTAL DATA RECONSTRUCTION SCRIPT (v1.6 - Discretization Fix)
%
% Description:
% This is the definitive script for reconstructing real experimental data.
% It closes the "sim-to-reality" gap by simulating the exact integer
% discretization of the signal profiles that occurs on the Arduino.
%
% v1.6 Improvements:
% - H-matrix generation now uses the exact, discretized profile values.
% - All previous bug fixes and correct pin maps are included.
% =========================================================================
clear; clc; close all;
%% ===== 1. CONFIGURATION AND FILE SELECTION =====
fprintf('=== Experimental Data Reconstruction Script (v1.6) ===\n');

% --- Output Setup ---
timestamp = datestr(now, 'mmddyy_HHMMSS');
output_folder = fullfile('reconstruction_output_real_data', timestamp);
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('Saving all results to: %s/\n', output_folder);

% --- HARDWARE CONFIGURATION (Must match experiment) ---
params.disabled_pmuts_grid = [10, 11]; 
params.fixed_rx_indices_grid = [3, 14]; 
params.num_active_tx = 8;
non_rx_grid_pos = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16];
hardware_pins   = 26:39; 
params.GRID_TO_PIN_MAP = containers.Map(non_rx_grid_pos, hardware_pins);

% --- ADMM HYPERPARAMETER TUNING ---
params.lambda_tv_reg = 25.0;
params.rho_admm = 6.73;

% --- Core Physical and Simulation Parameters ---
params.c = 343;
params.fs = 1e6;
params.pmut_width_m = 0.002;
params.kerf_m = 0.005;
params.excitation_amplitude = 1e12;
params.arduino_clock_hz = 16e6; % For discretization simulation

% --- Imaging Grid Parameters ---
params.grid_x_width_m = 0.150;
params.grid_y_width_m = 0.150;
params.target_height_m = 0.150;
params.grid_step_m = 0.004;

% --- Pre-processing & Gating Parameters ---
params.filter_cutoff_hz = 70000;
params.filter_order = 4;
params.gate_start_time_s = (params.target_height_m / params.c) * 0.8;
params.gate_duration_s = 3.0e-3;

% --- ADMM Solver Parameters ---
params.admm_tol = 1.2e-5;
params.admm_max_iter = 50;
params.pcg_max_iter = 30;
params.pcg_tol = 1e-8;
params.assembly_chunk_size = 25;

%% ===== 2. LOAD EXPERIMENTAL DATA AND PROFILES =====
fprintf('\n--- Step 1: Loading Experimental Data and Profiles ---\n');
[h5_file, h5_path] = uigetfile('*.h5', 'Select the TARGET HDF5 data file');
if isequal(h5_file, 0), error('User canceled file selection.'); end
data_filepath = fullfile(h5_path, h5_file);

[bg_h5_file, bg_h5_path] = uigetfile('*.h5', 'Select the corresponding BACKGROUND HDF5 file');
if isequal(bg_h5_file, 0), error('Background file is required for reconstruction.'); end
background_filepath = fullfile(bg_h5_path, bg_h5_file);

fprintf('  Loading data from HDF5 files...\n');
target_ch_a_raw = h5read(data_filepath, '/echo_data_ch_A');
target_ch_b_raw = h5read(data_filepath, '/echo_data_ch_B');
background_ch_a_raw = h5read(background_filepath, '/echo_data_ch_A');
background_ch_b_raw = h5read(background_filepath, '/echo_data_ch_B');
tx_pin_profiles = h5read(data_filepath, '/tx_pin_profiles'); % These are the HARDWARE PINS used

fprintf('  Please select the 4 corresponding CSV profile files...\n');
[delays_file, prof_path] = uigetfile('*.csv', 'Select DELAYS_US CSV file');
profiles.delays = readmatrix(fullfile(prof_path, delays_file));
[freqs_file, ~] = uigetfile('*.csv', 'Select FREQUENCIES_HZ CSV file', prof_path);
profiles.frequencies = readmatrix(fullfile(prof_path, freqs_file));
[phases_file, ~] = uigetfile('*.csv', 'Select PHASES_RAD CSV file', prof_path);
profiles.phases = readmatrix(fullfile(prof_path, phases_file));
[apods_file, ~] = uigetfile('*.csv', 'Select APODIZATIONS CSV file', prof_path);
profiles.apodizations = readmatrix(fullfile(prof_path, apods_file));

params.num_acquisitions = size(tx_pin_profiles, 2);
fprintf('  Loaded data for %d acquisitions.\n', params.num_acquisitions);

%% ===== 3. PRE-PROCESS EXPERIMENTAL DATA =====
% (This section is unchanged)
fprintf('\n--- Step 2: Pre-processing Experimental Data ---\n');
VOLTAGE_RANGE_MV = 500.0;
RESOLUTION_BITS = 14;
max_adc_value = 2^(RESOLUTION_BITS - 1) - 1;
target_a_mv = (double(target_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
target_b_mv = (double(target_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_a_mv = (double(background_ch_a_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
background_b_mv = (double(background_ch_b_raw) / max_adc_value) * VOLTAGE_RANGE_MV;
try
    timebase = h5readatt(data_filepath, '/echo_data_ch_A', 'timebase');
catch
    fprintf('  WARNING: "timebase" attribute not found. Defaulting to 4.\n');
    timebase = 4;
end
PICO_CLOCK_FREQ_HZ = 62.5e6;
if timebase < 3, sample_interval_s = (2^timebase) / (PICO_CLOCK_FREQ_HZ * 2);
else, sample_interval_s = (timebase - 2) / PICO_CLOCK_FREQ_HZ; end
original_fs = 1 / sample_interval_s;
fprintf('  Applying high-pass filter...\n');
[b, a] = butter(params.filter_order, params.filter_cutoff_hz / (original_fs / 2), 'high');
filtered_target_a = filtfilt(b, a, target_a_mv);
filtered_target_b = filtfilt(b, a, target_b_mv);
filtered_background_a = filtfilt(b, a, background_a_mv);
filtered_background_b = filtfilt(b, a, background_b_mv);
fprintf('  Performing background subtraction...\n');
subtracted_a = filtered_target_a - filtered_background_a;
subtracted_b = filtered_target_b - filtered_background_b;
final_a = subtracted_a - mean(subtracted_a, 1);
final_b = subtracted_b - mean(subtracted_b, 1);
fprintf('  Averaging receiver channels...\n');
processed_data = (final_a + final_b) / 2.0;
fprintf('  Resampling data to %.1f MHz...\n', params.fs / 1e6);
[p, q] = rat(params.fs / original_fs);
b_vector_resampled = resample(processed_data, p, q);
b_vector = b_vector_resampled(:);

%% ===== 4. GENERATE DIGITAL TWIN H-MATRIX =====
fprintf('\n--- Step 3: Generating Digital Twin H-Matrix ---\n');
tic;
[H_raw, imaging_grid] = generate_h_matrix(params, profiles, tx_pin_profiles);
fprintf('H-Matrix generation complete. Time: %.2f seconds.\n', toc);

% (The rest of the script: alignment, coherence, reconstruction is unchanged)
% ...
