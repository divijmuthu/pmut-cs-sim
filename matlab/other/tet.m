% --- Diagnostic Snippet ---
[file, path] = uigetfile('*.h5', 'Select your BACKGROUND HDF5 file to inspect');
if file
    filepath = fullfile(path, file);
    fprintf('--- Contents of %s ---\n', file);
    h5disp(filepath);
end