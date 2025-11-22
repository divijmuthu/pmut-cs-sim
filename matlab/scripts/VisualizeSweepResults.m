% VisualizeSweepResults.m
% Visualizes results from TheoryTest_Diagnostics.m parameter sweeps

clearvars; clc; close all;

results_dir = 'results';
visual_dir = 'results_visualized';
if ~exist(visual_dir, 'dir')
    mkdir(visual_dir);
end
files = dir(fullfile(results_dir, '*.mat'));

if isempty(files)
    error('No results found in %s', results_dir);
end

% Try to infer which parameter is being swept from filenames
param_names = {'lambda_tv', 'numItersADMM', 'R_acquisitions', 'target_SNR_db'};

for p = 1:numel(param_names)
    param = param_names{p};
    param_files = files(contains({files.name}, param));
    if isempty(param_files), continue; end

    param_vals = zeros(1, numel(param_files));
    final_psnr = zeros(1, numel(param_files));
    final_residual = zeros(1, numel(param_files));

    fig1 = figure('Name', ['Sweep: ' param], 'Position', [100 100 1200 600]);
    for k = 1:numel(param_files)
        data = load(fullfile(results_dir, param_files(k).name));
        % Extract parameter value from filename
        tokens = regexp(param_files(k).name, '([\-\d\.]+)', 'tokens');
        param_vals(k) = str2double(tokens{end}{1});
        final_psnr(k) = data.psnr(end);
        final_residual(k) = data.residuals(end);

        % Plot reconstruction
        subplot(2, numel(param_files), k);
        % --- Begin reshape logic ---
        if isfield(data, 'params') && isfield(data.params, 'grid_width_mm') && isfield(data.params, 'grid_step_mm')
            grid_width = data.params.grid_width_mm;
            grid_step = data.params.grid_step_mm;
            x_coords = -grid_width/2 : grid_step : grid_width/2;
            y_coords = -grid_width/2 : grid_step : grid_width/2;
            img_size = [length(y_coords), length(x_coords)];
            if numel(data.reconstruction) == prod(img_size)
                recon_img = reshape(data.reconstruction, img_size);
                imagesc(x_coords, y_coords, recon_img);
                xlabel('x (mm)'); ylabel('y (mm)');
            else
                imagesc(data.reconstruction); % fallback
            end
        else
            % fallback: try to guess square
            N = sqrt(length(data.reconstruction));
            if mod(N,1)==0
                imagesc(reshape(data.reconstruction, [N N]));
            else
                imagesc(data.reconstruction); % fallback
            end
        end
        % --- End reshape logic ---
        axis image off; colormap gray;
        title(sprintf('%s = %.2g', param, param_vals(k)));
        if isfield(data, 'scene_matrix')
            hold on; contour(data.scene_matrix, [0.5 0.5], 'r', 'LineWidth', 1); hold off;
        end

        % Plot PSNR/residual curves
        subplot(2, numel(param_files), k+numel(param_files));
        yyaxis left; plot(data.psnr, 'b-o'); ylabel('PSNR (dB)');
        yyaxis right; plot(log10(data.residuals), 'r-x'); ylabel('log_{10}(Residual)');
        xlabel('Iteration');
        title('Convergence');
        grid on;

        % Save each reconstruction+convergence plot
        saveas(fig1, fullfile(visual_dir, sprintf('Sweep_%s_%g.png', param, param_vals(k))));
    end

    % Summary plot
    fig2 = figure('Name', ['Summary: ' param], 'Position', [400 100 600 400]);
    subplot(2,1,1);
    plot(param_vals, final_psnr, 'bo-'); xlabel(param); ylabel('Final PSNR (dB)'); grid on;
    subplot(2,1,2);
    plot(param_vals, final_residual, 'rx-'); xlabel(param); ylabel('Final Residual'); grid on;
    saveas(fig2, fullfile(visual_dir, sprintf('Summary_%s.png', param)));
end

disp(['Visualization complete. Plots saved in ' visual_dir]); 