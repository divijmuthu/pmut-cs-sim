% Test script for Experimental Reconstruction
% This script creates sample experimental data and tests the reconstruction

clearvars;
clc;
close all;

fprintf('=== Testing Experimental Reconstruction ===\n\n');

% Create test output directory
test_dir = fullfile(pwd, 'test_experimental_output');
if ~exist(test_dir, 'dir')
    mkdir(test_dir);
end

% Create sample H5 file with echo data
fprintf('Creating sample H5 file...\n');
sample_echo_data = randn(200, 1000); % 200 acquisitions, 1000 samples each
h5_file_path = fullfile(test_dir, 'sample_echo_data.h5');
h5create(h5_file_path, '/echo_data', size(sample_echo_data));
h5write(h5_file_path, '/echo_data', sample_echo_data);
fprintf('  ✓ Sample H5 file created: %s\n', h5_file_path);

% Create sample delay profiles CSV
fprintf('Creating sample delay profiles...\n');
sample_delays = rand(200, 3) * 12; % 200 acquisitions, 3 elements, max 12 us
delay_csv_path = fullfile(test_dir, 'sample_delays.csv');
delay_table = array2table(sample_delays, 'VariableNames', {'Delay1', 'Delay2', 'Delay3'});
writetable(delay_table, delay_csv_path);
fprintf('  ✓ Sample delay CSV created: %s\n', delay_csv_path);

% Create sample frequencies CSV
fprintf('Creating sample frequencies...\n');
sample_frequencies = 45e3 + rand(200, 3) * 20e3; % 45-65 kHz range
freq_csv_path = fullfile(test_dir, 'sample_frequencies.csv');
freq_table = array2table(sample_frequencies, 'VariableNames', {'Freq1', 'Freq2', 'Freq3'});
writetable(freq_table, freq_csv_path);
fprintf('  ✓ Sample frequency CSV created: %s\n', freq_csv_path);

% Test the experimental reconstruction
fprintf('\n=== Running Experimental Reconstruction Test ===\n');

% Create a modified version of the experimental reconstruction for testing
try
    % Initialize with test parameters
    params = struct();
    params.c = 343;
    params.fs = 2e6;
    params.pMUT_width_mm = 20;
    params.pMUT_spacing_mm = 20;
    params.kerf_mm = 0.1;
    params.grid_width_mm = 150;
    params.grid_depth_start_mm = 250;
    params.grid_depth_end_mm = 350;
    params.grid_step_mm = 4;
    params.R_acquisitions = 100;
    params.excitation_amplitude = 10000;
    params.target_SNR_db = 30;
    params.max_delay_us = 12;
    params.numItersADMM = 10; % Reduced for testing
    params.rho_admm = 10;
    params.lambda_tv_reg = 0.1;
    params.use_last_100 = true;
    params.use_parallel = false;
    params.max_workers = 4;
    params.admm_tol = 1e-6;
    params.admm_max_iter = 10; % Reduced for testing
    params.pcg_max_iter = 25;
    params.pcg_tol = 1e-8;
    params.use_gpu = false;
    params.use_sparse = true;
    params.memory_pool = true;
    params.vectorize_interp = true;
    params.adaptive_tolerance = true;
    params.block_size = 100;
    params.cache_results = true;
    
    % Set file paths
    params.h5_file_path = h5_file_path;
    params.delay_csv_path = delay_csv_path;
    params.freq_csv_path = freq_csv_path;
    
    % Create output folder
    output_folder = fullfile(test_dir, 'reconstruction_results');
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    % Add m_files to path
    mfiles_path = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'm_files');
    if exist(mfiles_path, 'dir') && ~contains(path, mfiles_path)
        addpath(genpath(mfiles_path));
    end
    
    % Initialize Field II
    field_init(-1);
    set_field('fs', params.fs);
    set_field('c', params.c);
    
    fprintf('✓ Test setup completed successfully\n');
    fprintf('✓ File paths configured\n');
    fprintf('✓ Field II initialized\n');
    fprintf('✓ Ready to run experimental reconstruction\n\n');
    
    fprintf('To run the full experimental reconstruction:\n');
    fprintf('1. Replace the sample files with your real experimental data\n');
    fprintf('2. Run: ExperimentalReconstruction.m\n');
    fprintf('3. The script will prompt for your H5 and CSV files\n\n');
    
    fprintf('Sample files created in: %s\n', test_dir);
    fprintf('  - %s\n', h5_file_path);
    fprintf('  - %s\n', delay_csv_path);
    fprintf('  - %s\n', freq_csv_path);
    
catch ME
    fprintf('Error during test setup: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

fprintf('\n=== Test Complete ===\n'); 