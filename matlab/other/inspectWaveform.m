% =========================================================================
% MAT File Inspector Script
%
% Description:
% This script helps diagnose issues with PicoScope .mat files. It prompts
% the user to select a single .mat file, then loads it and displays the
% names, sizes, and classes of all variables stored inside. This helps
% identify the correct variable name and orientation of the waveform data.
% =========================================================================
clear; clc; close all;

fprintf('=== PicoScope .mat File Inspector ===\n\n');

% 1. Select a single .mat file from the exported folder
[mat_file, mat_path] = uigetfile('*.mat', 'Select a SINGLE .mat file from your waveform folder');
if isequal(mat_file, 0)
    disp('Operation canceled by user.');
    return;
end
filepath = fullfile(mat_path, mat_file);
fprintf('Inspecting file: %s\n\n', filepath);

% 2. Load the data into a struct
try
    data_struct = load(filepath);
catch ME
    fprintf('Error loading file: %s\n', ME.message);
    return;
end

% 3. Analyze and display the contents
fprintf('--- Contents of the .mat file ---\n');
variable_names = fieldnames(data_struct);

for i = 1:length(variable_names)
    var_name = variable_names{i};
    var_content = data_struct.(var_name);
    var_class = class(var_content);
    var_size = size(var_content);
    
    fprintf('Variable Name: "%s"\n', var_name);
    fprintf('  - Class: %s\n', var_class);
    fprintf('  - Size: %s\n', mat2str(var_size));
    
    % Heuristic to identify the likely waveform data
    if isnumeric(var_content) && isvector(var_content) && length(var_content) > 1000
        fprintf('  - LIKELY WAVEFORM DATA\n');
    end
    fprintf('\n');
end

fprintf('--- Analysis Complete ---\n');
fprintf('Please look for the variable marked as "LIKELY WAVEFORM DATA".\n');
fprintf('Note its name and size (e.g., 1x125004 for a row vector).\n');
