%% analyze_v29_actual_results.m - Analyze Actual V29_realistic Results
% Compare actual results with V28_final and analyze realistic pMUT performance

clear; clc; close all;

%% Load actual V29_realistic results
csv_file = 'sweep_output_realistic_v29/072525_165125/realistic_summary.csv';
fprintf('=== Analyzing Actual V29_realistic Results ===\n');
fprintf('Loading data from: %s\n\n', csv_file);

data = readtable(csv_file);
num_runs = height(data);

fprintf('Loaded %d completed runs with actual metrics\n\n', num_runs);

%% Display actual results summary
fprintf('=== ACTUAL V29_realistic Results Summary ===\n');

% Extract key metrics
max_coherence = data.max_coherence;
condition_numbers = data.condition_number;
rip_proxies = data.rip_proxy;
sparsity = data.sparsity;

fprintf('Performance Ranges:\n');
fprintf('  Max Coherence: %.3f - %.3f\n', min(max_coherence), max(max_coherence));
fprintf('  Condition Number: %.1f - %.1f\n', min(condition_numbers), max(condition_numbers));
fprintf('  RIP Proxy: %.3f - %.3f\n', min(rip_proxies), max(rip_proxies));
fprintf('  Sparsity: %.1f%% - %.1f%%\n', min(sparsity), max(sparsity));

%% Performance assessment
fprintf('\n=== Performance Assessment ===\n');

% Coherence assessment
excellent_coherence = sum(max_coherence < 0.85);
good_coherence = sum(max_coherence < 0.9);
poor_coherence = sum(max_coherence >= 0.9);

fprintf('Coherence Performance:\n');
fprintf('  🏆 OUTSTANDING (<0.85): %d/%d runs (%.1f%%)\n', excellent_coherence, num_runs, excellent_coherence/num_runs*100);
fprintf('  ✅ EXCELLENT (<0.9): %d/%d runs (%.1f%%)\n', good_coherence, num_runs, good_coherence/num_runs*100);
fprintf('  ⚠️ POOR (>=0.9): %d/%d runs (%.1f%%)\n', poor_coherence, num_runs, poor_coherence/num_runs*100);

% Condition number assessment
excellent_condition = sum(condition_numbers < 100);
good_condition = sum(condition_numbers < 200);
poor_condition = sum(condition_numbers >= 200);

fprintf('\nCondition Number Performance:\n');
fprintf('  🏆 OUTSTANDING (<100): %d/%d runs (%.1f%%)\n', excellent_condition, num_runs, excellent_condition/num_runs*100);
fprintf('  ✅ EXCELLENT (<200): %d/%d runs (%.1f%%)\n', good_condition, num_runs, good_condition/num_runs*100);
fprintf('  ⚠️ POOR (>=200): %d/%d runs (%.1f%%)\n', poor_condition, num_runs, poor_condition/num_runs*100);

%% Best performers analysis
fprintf('\n=== BEST PERFORMERS (Actual Results) ===\n');

% Find best coherence
[best_coherence, best_coherence_idx] = min(max_coherence);
[best_condition, best_condition_idx] = min(condition_numbers);
[best_rip, best_rip_idx] = min(rip_proxies);

fprintf('Best Coherence: Run %d (%.4f)\n', data.run_index(best_coherence_idx), best_coherence);
fprintf('  Parameters: Impulse=%dμs, Foff=%dHz, Acq=%d\n', ...
    data.impulse_duration_us(best_coherence_idx), data.frequency_offset_hz(best_coherence_idx), ...
    data.num_acquisitions(best_coherence_idx));

fprintf('Best Condition Number: Run %d (%.2f)\n', data.run_index(best_condition_idx), best_condition);
fprintf('  Parameters: Impulse=%dμs, Foff=%dHz, Acq=%d\n', ...
    data.impulse_duration_us(best_condition_idx), data.frequency_offset_hz(best_condition_idx), ...
    data.num_acquisitions(best_condition_idx));

fprintf('Best RIP Proxy: Run %d (%.4f)\n', data.run_index(best_rip_idx), best_rip);
fprintf('  Parameters: Impulse=%dμs, Foff=%dHz, Acq=%d\n', ...
    data.impulse_duration_us(best_rip_idx), data.frequency_offset_hz(best_rip_idx), ...
    data.num_acquisitions(best_rip_idx));

%% Parameter impact analysis
fprintf('\n=== Parameter Impact Analysis ===\n');

% Impulse duration impact
unique_imp = unique(data.impulse_duration_us);
fprintf('Impulse Duration Impact:\n');
for imp = unique_imp'
    runs_with_imp = data.impulse_duration_us == imp;
    avg_coherence = mean(data.max_coherence(runs_with_imp));
    avg_condition = mean(data.condition_number(runs_with_imp));
    fprintf('  %dμs: Coherence=%.3f, Condition=%.1f (%d runs)\n', ...
        imp, avg_coherence, avg_condition, sum(runs_with_imp));
end

% Frequency offset impact
unique_foff = unique(data.frequency_offset_hz);
fprintf('\nFrequency Offset Impact:\n');
for foff = unique_foff'
    runs_with_foff = data.frequency_offset_hz == foff;
    avg_coherence = mean(data.max_coherence(runs_with_foff));
    avg_condition = mean(data.condition_number(runs_with_foff));
    fprintf('  %dHz: Coherence=%.3f, Condition=%.1f (%d runs)\n', ...
        foff, avg_coherence, avg_condition, sum(runs_with_foff));
end

% Acquisition count impact
unique_acq = unique(data.num_acquisitions);
fprintf('\nAcquisition Count Impact:\n');
for acq = unique_acq'
    runs_with_acq = data.num_acquisitions == acq;
    avg_coherence = mean(data.max_coherence(runs_with_acq));
    avg_condition = mean(data.condition_number(runs_with_acq));
    fprintf('  %d acq: Coherence=%.3f, Condition=%.1f (%d runs)\n', ...
        acq, avg_coherence, avg_condition, sum(runs_with_acq));
end

%% Comparison with V28_final
fprintf('\n=== Comparison with V28_final (Idealized Waveforms) ===\n');

% V28_final best results
v28_hanning_coherence = [0.5647, 0.5796, 0.5788];
v28_hanning_condition = [20.51, 18.22, 18.44];
v28_gaussian_coherence = [0.8893, 0.8684, 0.8513];
v28_gaussian_condition = [84.33, 84.76, 80.74];
v28_chirp_coherence = [0.9061, 0.8783, 0.8768];
v28_chirp_condition = [95.55, 77.59, 85.50];

fprintf('V28_final Best Results:\n');
fprintf('  Hanning pattern: Coherence=%.3f-%.3f, Condition=%.1f-%.1f\n', ...
    min(v28_hanning_coherence), max(v28_hanning_coherence), ...
    min(v28_hanning_condition), max(v28_hanning_condition));
fprintf('  Gaussian pattern: Coherence=%.3f-%.3f, Condition=%.1f-%.1f\n', ...
    min(v28_gaussian_coherence), max(v28_gaussian_coherence), ...
    min(v28_gaussian_condition), max(v28_gaussian_condition));
fprintf('  Chirp pattern: Coherence=%.3f-%.3f, Condition=%.1f-%.1f\n', ...
    min(v28_chirp_coherence), max(v28_chirp_coherence), ...
    min(v28_chirp_condition), max(v28_chirp_condition));

fprintf('\nV29_realistic Actual Results:\n');
fprintf('  Realistic impulse: Coherence=%.3f-%.3f, Condition=%.1f-%.1f\n', ...
    min(max_coherence), max(max_coherence), min(condition_numbers), max(condition_numbers));

% Calculate degradation
v28_best_coherence = min(v28_hanning_coherence);
v29_best_coherence = min(max_coherence);
coherence_degradation = (v29_best_coherence - v28_best_coherence) / v28_best_coherence * 100;

v28_best_condition = min(v28_hanning_condition);
v29_best_condition = min(condition_numbers);
condition_degradation = (v29_best_condition - v28_best_condition) / v28_best_condition * 100;

fprintf('\nPerformance Impact of Realistic Constraints:\n');
fprintf('  Coherence degradation: %.1f%% (%.3f → %.3f)\n', ...
    coherence_degradation, v28_best_coherence, v29_best_coherence);
fprintf('  Condition number degradation: %.1f%% (%.1f → %.1f)\n', ...
    condition_degradation, v28_best_condition, v29_best_condition);

%% Optimal realistic parameters
fprintf('\n=== Optimal Realistic Parameters ===\n');

% Find best overall performer (combined score)
combined_scores = max_coherence + condition_numbers/1000;
[best_score, best_score_idx] = min(combined_scores);

best_run = data(best_score_idx, :);
fprintf('Best Overall Performer: Run %d (Score=%.4f)\n', best_run.run_index, best_score);
fprintf('  Max Coherence: %.4f\n', best_run.max_coherence);
fprintf('  Condition Number: %.2f\n', best_run.condition_number);
fprintf('  RIP Proxy: %.4f\n', best_run.rip_proxy);
fprintf('  Parameters: Impulse=%dμs, Foff=%dHz, Acq=%d\n', ...
    best_run.impulse_duration_us, best_run.frequency_offset_hz, best_run.num_acquisitions);

%% Detailed results table
fprintf('\n=== Detailed Results Table ===\n');
fprintf('Run | Impulse | Foff | Acq | Coherence | Condition | RIP | Performance\n');
fprintf('----|---------|------|-----|-----------|-----------|-----|------------\n');

for i = 1:num_runs
    run_data = data(i, :);
    
    % Performance assessment
    if run_data.max_coherence < 0.85
        perf = '🏆';
    elseif run_data.max_coherence < 0.9
        perf = '✅';
    else
        perf = '⚠️';
    end
    
    fprintf('%3d | %7dμs | %4d | %3d | %9.4f | %9.2f | %.3f | %s\n', ...
        run_data.run_index, run_data.impulse_duration_us, run_data.frequency_offset_hz, ...
        run_data.num_acquisitions, run_data.max_coherence, run_data.condition_number, ...
        run_data.rip_proxy, perf);
end

%% Key insights
fprintf('\n=== Key Insights ===\n');

fprintf('1. Realistic pMUT Performance:\n');
fprintf('   - Best coherence: %.4f (vs %.4f with idealized waveforms)\n', v29_best_coherence, v28_best_coherence);
fprintf('   - Best condition: %.1f (vs %.1f with idealized waveforms)\n', v29_best_condition, v28_best_condition);
fprintf('   - Performance degradation: %.1f%% coherence, %.1f%% condition\n', coherence_degradation, condition_degradation);

fprintf('\n2. Parameter Optimization:\n');
best_imp = data.impulse_duration_us(best_coherence_idx);
best_foff = data.frequency_offset_hz(best_coherence_idx);
best_acq = data.num_acquisitions(best_coherence_idx);
fprintf('   - Optimal impulse duration: %dμs\n', best_imp);
fprintf('   - Optimal frequency offset: %dHz\n', best_foff);
fprintf('   - Optimal acquisition count: %d\n', best_acq);

fprintf('\n3. Implementation Strategy:\n');
fprintf('   - Use short impulse excitation (%dμs) at 57.7 kHz\n', best_imp);
fprintf('   - Maintain time-domain diversity (interpolation approach)\n');
fprintf('   - Accept %.1f%% performance degradation for realistic constraints\n', coherence_degradation);

%% Save comprehensive analysis
fprintf('\n=== Saving Comprehensive Analysis ===\n');

analysis_summary = struct();
analysis_summary.v29_best_coherence = v29_best_coherence;
analysis_summary.v29_best_condition = v29_best_condition;
analysis_summary.v28_best_coherence = v28_best_coherence;
analysis_summary.v28_best_condition = v28_best_condition;
analysis_summary.coherence_degradation = coherence_degradation;
analysis_summary.condition_degradation = condition_degradation;
analysis_summary.optimal_parameters = struct('impulse_duration', best_imp, 'frequency_offset', best_foff, 'acquisition_count', best_acq);
analysis_summary.all_data = data;

save(fullfile('sweep_output_realistic_v29/072525_165125', 'v29_comprehensive_analysis.mat'), 'analysis_summary');

fprintf('Comprehensive analysis saved to: v29_comprehensive_analysis.mat\n');

fprintf('\n=== V29_realistic Actual Results Analysis Complete! ===\n');
fprintf('Key Finding: Realistic pMUT constraints cause %.1f%% coherence degradation\n', coherence_degradation);
fprintf('but still achieve excellent performance with proper parameter optimization.\n'); 