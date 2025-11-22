% Quick H Matrix Check
% Check if the H matrix is now non-zero after the time delay fix
clearvars;
clc;

%% Load the latest H matrix
fprintf('--- Checking Latest H Matrix ---\n');

h_matrix_file = '/Users/deepshikhakaul/Documents/MATLAB/Optimized_071925_012/optimized_H_matrix.mat';
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
        
        % Check condition number
        cond_num = cond(H);
        fprintf('  Condition number: %.2e\n', cond_num);
        
        if sum(H(:) ~= 0) > 0
            fprintf('✅ H matrix is now non-zero!\n');
        else
            fprintf('❌ H matrix is still all zeros!\n');
        end
    else
        error('No H matrix found in file');
    end
else
    error('H matrix file not found');
end

%% Check a few sample rows
fprintf('\n--- Sample Row Analysis ---\n');
for i = 1:min(5, size(H, 1))
    row = H(i, :);
    fprintf('Row %d: mean=%.6f, std=%.6f, min=%.6f, max=%.6f, non-zero=%d\n', ...
        i, mean(row), std(row), min(row), max(row), sum(row ~= 0));
end

fprintf('\n--- Check Complete ---\n'); 