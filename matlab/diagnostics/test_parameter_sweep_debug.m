%% test_parameter_sweep_debug.m - Debug parameter sweep issues
% Test individual components to find the problem

clear; clc; close all;

fprintf('=== DEBUGGING PARAMETER SWEEP ===\n');

% Test basic configuration
config = struct();
config.c = 343;
config.fs = 1e6;
config.pmut_width_m = 0.020;
config.tx_pool_width_m = 0.200;
config.grid_width_m = 0.150;
config.target_distance_m = 0.150;
config.grid_depth_range_m = 0.100;
config.grid_step_m = 0.0075;  % 7.5mm
config.excitation_amplitude = 1e15;
config.pmut_resonance_freq = 57700;
config.pmut_bandwidth = 2520;
config.impulse_duration_us = 10;
config.num_active_tx = 5;
config.max_delay_rand_us = 500;
config.apodization_mode = 'uniform';
config.frequency_offset_hz = 0;
config.num_acquisitions = 20;
config.numItersADMM = 50;
config.rho_admm = 6.73;
config.lambda_tv_reg = 0.7438;
config.admm_tol = 1.2e-5;
config.admm_max_iter = 50;
config.pcg_max_iter = 30;
config.pcg_tol = 1e-8;
config.target_SNR_db = 35;

% Test parameters
params = struct();
params.grid_step_mm = 7.5;
params.target_size_mm = 4;
params.grid_spacing_mm = 20;
params.num_acquisitions = 20;
params.lambda_tv_reg = 0.7438;
params.target_pattern = '3x3_grid';

fprintf('Testing H matrix generation...\n');
try
    [H, coherence_matrix] = generate_h_matrix_realistic_pmut(config);
    fprintf('✓ H matrix generation successful\n');
    fprintf('  H matrix size: %d x %d\n', size(H));
    fprintf('  Non-zero elements: %d\n', nnz(H));
catch ME
    fprintf('✗ H matrix generation failed: %s\n', ME.message);
    return;
end

fprintf('\nTesting scene creation...\n');
try
    % Create imaging grid
    x_coords_img = -config.grid_width_m/2 : config.grid_step_m : config.grid_width_m/2;
    z_coords_img = (config.target_distance_m - config.grid_depth_range_m/2) : config.grid_step_m : (config.target_distance_m + config.grid_depth_range_m/2);
    [X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
    
    scene_matrix = create_scene_with_pattern(X_mesh, Z_mesh, config, params);
    fprintf('✓ Scene creation successful\n');
    fprintf('  Scene size: %d x %d\n', size(scene_matrix));
    fprintf('  Non-zero pixels: %d\n', nnz(scene_matrix));
catch ME
    fprintf('✗ Scene creation failed: %s\n', ME.message);
    return;
end

fprintf('\nTesting forward simulation...\n');
try
    v_true_vector = scene_matrix(:);
    Hv_signal = H * v_true_vector;
    signal_power = mean(Hv_signal(:).^2);
    target_SNR_linear = 10^(config.target_SNR_db / 10);
    noise_variance = signal_power / target_SNR_linear;
    noise_sigma = sqrt(noise_variance);
    noise = noise_sigma * randn(size(Hv_signal));
    b_vector = Hv_signal + noise;
    fprintf('✓ Forward simulation successful\n');
    fprintf('  Signal power: %.2e\n', signal_power);
    fprintf('  Noise variance: %.2e\n', noise_variance);
catch ME
    fprintf('✗ Forward simulation failed: %s\n', ME.message);
    return;
end

fprintf('\nTesting ADMM reconstruction...\n');
try
    output_folder = 'debug_test_output';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    reconstructed_image = run_admm_reconstruction(H, b_vector, scene_matrix, config, output_folder);
    fprintf('✓ ADMM reconstruction successful\n');
    fprintf('  Reconstruction size: %d x %d\n', size(reconstructed_image));
catch ME
    fprintf('✗ ADMM reconstruction failed: %s\n', ME.message);
    fprintf('  Error details: %s\n', getReport(ME, 'extended'));
    return;
end

fprintf('\nTesting metrics calculation...\n');
try
    scene_norm = scene_matrix / max(scene_matrix(:));
    recon_norm = reconstructed_image / max(reconstructed_image(:));
    
    MSE = mean((scene_norm(:) - recon_norm(:)).^2);
    psnr = 10 * log10(1 / MSE);
    
    scene_vec = scene_norm(:);
    recon_vec = recon_norm(:);
    correlation = (scene_vec - mean(scene_vec))' * (recon_vec - mean(recon_vec)) / ...
        (sqrt(sum((scene_vec - mean(scene_vec)).^2)) * sqrt(sum((recon_vec - mean(recon_vec)).^2)));
    
    fprintf('✓ Metrics calculation successful\n');
    fprintf('  PSNR: %.2f dB\n', psnr);
    fprintf('  Correlation: %.4f\n', correlation);
catch ME
    fprintf('✗ Metrics calculation failed: %s\n', ME.message);
    return;
end

fprintf('\n=== DEBUG TEST COMPLETE ===\n');
fprintf('All components working correctly!\n'); 