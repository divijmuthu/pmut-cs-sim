% Advanced_pMUT_CS_Recon_v6.0.0_TimeDelayCoding.m
%
% Description:
% This version marks a significant change in the H-matrix generation strategy,
% moving from phase-coding to time-delay coding.
%
% Key Changes:
% 1. The high-energy chirp burst is now defined as the transducer's primary
%    impulse response via xdc_impulse().
% 2. The acquisition loop now uses xdc_focus_times() to apply a random
%    time delay to each of the 3 active pMUTs for each acquisition.
% 3. The ele_waveform() method is no longer used for coding.
%
% This approach is inspired by the successful coding mechanism in the
% reference Ensemble.m script and aims to improve H-matrix conditioning.
% All other parameters (geometry, grid, SNR, reconstruction algorithms)
% are kept consistent with v5.7.0 for direct comparison.
%
clearvars; clc; close all;

% --- Initialize Field II ---
field_init(-1);

% --- Core Physical and Simulation Constants ---
c = 1540;           % Speed of Sound [m/s]
fs = 2e6;           % Simulation Sampling Frequency [Hz].
fc_nominal = 1.0e5; % Reference nominal Center Frequency for chirp (100 kHz)
lambda = c/fc_nominal;
set_field('fs', fs);
set_field('c', c);

fprintf('--- v6.0.0: Time-Delay Coding for H-Matrix Enhancement ---\n');
fprintf('\n');

% --- pMUT Element and Array Geometry ---
pMUT_width_mm = 20;         % Active pMUT width/height (mm)
pMUT_spacing_mm = 20;       % Triangle sides (mm)
kerf_mm = 0.1;              % Kerf between virtual grid elements

% --- Imaging Grid Geometry ---
grid_width_mm = 150;        % Lateral imaging width (mm)
grid_depth_start_mm = 250;  % Start imaging at (mm)
grid_depth_end_mm = 350;    % End imaging at (mm)
grid_step_mm = 4;           % Pixel size (mm)

% --- SIMULATION CONTROL PARAMETERS ---
R_acquisitions = 50; % more ofc slower but marked improvement!
excitation_amplitude = 500;
electronic_SNR_db = 90; % not bad even with 60 lol
max_delay_us = 20;          % Maximum random delay in microseconds for coding

% --- Reconstruction Algorithm Parameters (Defaults) ---
maxItersCG_main = 1000;
tolCG_main = 1e-8;
numItersADMM = 100;
rho_admm = 10;
lambda_tv_reg = 0.1;
% l_tv_reg = 0.1 produced excellent results
% more ADMM iters of 200 looks good but slow lol
% lowered rho_admm to 1, faster to reach but worse peak, 10 was fine, 50
% not much better eiether

fprintf('--- Key Simulation Parameters ---\n');
fprintf('pMUT Width: %g mm, Spacing: %g mm\n', pMUT_width_mm, pMUT_spacing_mm);
fprintf('Imaging Grid: %.0f-%.0f mm depth, %.0f mm width, %.1f mm step\n', grid_depth_start_mm, grid_depth_end_mm, grid_width_mm, grid_step_mm);
fprintf('Acquisitions: %d, Max Delay: %g us, SNR: %g dB\n', R_acquisitions, max_delay_us, electronic_SNR_db);
fprintf('\n');

% --- Convert mm parameters to meters for Field II ---
pMUT_width = pMUT_width_mm/1000;
pMUT_height = pMUT_width; % Assume square
kerf = kerf_mm/1000;
d_spacing = pMUT_spacing_mm/1000;

grid_width = grid_width_mm/1000;
grid_depth_start = grid_depth_start_mm/1000;
grid_depth_end = grid_depth_end_mm/1000;
grid_step = grid_step_mm/1000;
max_delay = max_delay_us / 1e6; % Convert max delay to seconds

% --- Define Imaging Grid (Voxels) ---
x_coords_img = -grid_width/2 : grid_step : grid_width/2;
z_coords_img = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
Y_mesh = zeros(size(X_mesh));
N_pixels = numel(X_mesh);
fprintf('Total imaging pixels (N_pixels): %d (%d axial x %d lateral).\n', N_pixels, length(z_coords_img), length(x_coords_img));
fprintf('\n');

% --- Define pMUT Aperture ---
% (Identical logic to v5.7.0)
triangle_side_length = d_spacing;
R_circ = triangle_side_length / sqrt(3);
pos1 = [R_circ, 0, 0];
pos2 = [R_circ*cos(2*pi/3), R_circ*sin(2*pi/3), 0];
pos3 = [R_circ*cos(4*pi/3), R_circ*sin(4*pi/3), 0];
desired_positions = [pos1; pos2; pos3];
num_active_intended = size(desired_positions, 1);
fprintf('Intending to use %d pMUTs.\n', num_active_intended);
num_x_grid = 7; num_y_grid = 7;
fprintf('Using a %dx%d virtual grid for mapping.\n', num_x_grid, num_y_grid);
physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
element_no_grid_map = 0;
center_offset_x = (num_x_grid - 1)/2 * (pMUT_width + kerf);
center_offset_y = (num_y_grid - 1)/2 * (pMUT_height + kerf);
for iy = 1:num_y_grid; y_pos_el = (iy-1)*(pMUT_height+kerf) - center_offset_y; for ix = 1:num_x_grid; x_pos_el = (ix-1)*(pMUT_width+kerf) - center_offset_x; element_no_grid_map = element_no_grid_map+1; physical_element_centers(element_no_grid_map,:) = [x_pos_el,y_pos_el,0]; end; end
active_indices_linear = zeros(num_active_intended, 1);
for i=1:num_active_intended; distances=sum((physical_element_centers-desired_positions(i,:)).^2,2); [~,min_idx]=min(distances); active_indices_linear(i)=min_idx; end
active_indices_linear = unique(active_indices_linear);
num_active = length(active_indices_linear);
if num_active<num_active_intended; warning('Could only map %d unique elements.',num_active); end
if num_active==0; error('No active elements found!'); end
fprintf('Successfully mapped %d unique active pMUTs.\n', num_active);
enabled_matrix=zeros(num_y_grid,num_x_grid); [row_indices,col_indices]=ind2sub([num_y_grid,num_x_grid],active_indices_linear);
for i=1:num_active; enabled_matrix(row_indices(i),col_indices(i))=1; end
pMUT_Aperture = xdc_2d_array(num_x_grid,num_y_grid,pMUT_width,pMUT_height,kerf,kerf,enabled_matrix,1,1,[0 0 100e-3]);

% --- Define High-Energy Chirp as the IMPULSE RESPONSE ---
f_start_chirp = 10e3;
f_end_chirp = 200e3;
burst_duration = 0.02e-3; % 20 us
fprintf('Synthetic Burst: %g-%g kHz, duration %g us, amplitude %g\n', f_start_chirp/1e3, f_end_chirp/1e3, burst_duration*1e3, excitation_amplitude);
t_burst_vec = 0 : 1/fs : burst_duration;
synth_burst_base = chirp(t_burst_vec, f_start_chirp, t_burst_vec(end), f_end_chirp, 'linear');
synth_burst_windowed = synth_burst_base .* tukeywin(length(t_burst_vec), 0.25)';
impulse_response_waveform = synth_burst_windowed * excitation_amplitude;

% Set this burst as the impulse response for the aperture
xdc_impulse(pMUT_Aperture, impulse_response_waveform);

% The excitation is now a simple, single trigger pulse.
xdc_excitation(pMUT_Aperture, 1);
fprintf('\n');

figure(2); clf;
plot(t_burst_vec*1e6, impulse_response_waveform);
title('Figure 2: Excitation Burst (now used as impulse response)');
xlabel('Time (us)'); ylabel('Amplitude'); grid on; set(gcf, 'Color', 'w');


%% Generate Coded Acquisitions and Measure H components
fprintf('\n--- Generating H Matrix Components (Time-Delay Coding) ---\n');
hydrophone_positions_img = [X_mesh(:), Y_mesh(:), Z_mesh(:)];

all_hhp_data = cell(R_acquisitions, 1);
all_start_times = zeros(R_acquisitions, 1);
all_K_values = zeros(R_acquisitions, 1);
rng('default'); % For reproducible random delays

for r_acq = 1:R_acquisitions
    fprintf('Calculating Acquisition %d/%d...\n', r_acq, R_acquisitions);
    if num_active == 0; error('Cannot proceed with 0 active elements.'); end
    
    % --- Generate a vector of random time delays ---
    % One delay for each active element, from 0 to max_delay
    delay_vector = rand(1, num_active) * max_delay;
    
    % --- Apply delays to the aperture elements ---
    xdc_focus_times(pMUT_Aperture, 0, delay_vector);
    
    % --- Calculate the combined response for this coded acquisition ---
    [hhp_r, start_time_r] = calc_hhp(pMUT_Aperture, pMUT_Aperture, hydrophone_positions_img);
    
    all_hhp_data{r_acq} = hhp_r;
    all_start_times(r_acq) = start_time_r;
    all_K_values(r_acq) = size(hhp_r, 1);
    
    if ~isempty(hhp_r)
        fprintf('    Max abs value in hhp_r: %g (K_r = %d)\n', max(abs(hhp_r(:))), size(hhp_r,1));
    else
        fprintf('    WARNING: hhp_r is empty for acquisition %d!\n', r_acq);
    end
end
fprintf('\n');

% The rest of the script is identical to v5.7.0: H-assembly, scene generation,
% noise, PCG and ADMM reconstruction, and plotting.

%% Assemble Full H Matrix with Interpolation and Alignment
fprintf('\n--- Assembling Full H Matrix ---\n');
all_end_times=zeros(R_acquisitions,1); 
for r_idx=1:R_acquisitions 
    if all_K_values(r_idx)>0 
        all_end_times(r_idx)=all_start_times(r_idx)+(all_K_values(r_idx)-1)/fs; 
    end 
end
min_global_start_time=min(all_start_times); 
max_global_end_time=max(all_end_times);
if isempty(min_global_start_time) || isempty(max_global_end_time) || min_global_start_time>=max_global_end_time
    min_global_start_time=0; max_K_val=max(all_K_values(all_K_values>0)); 
    if isempty(max_K_val) || max_K_val == 0; max_K_val=100; end 
    max_global_end_time=(max_K_val-1)/fs; 
    if max_global_end_time < min_global_start_time; max_global_end_time = min_global_start_time + (100-1)/fs; end
    warning('H assembly time window invalid, using fallback.');
end
t_common_axis=min_global_start_time:1/fs:max_global_end_time; K_global_per_acq=length(t_common_axis);
if K_global_per_acq == 0; K_global_per_acq = 1; t_common_axis=min_global_start_time; warning('K_global_per_acq was 0, set to 1.'); end
fprintf('Global Time Window for H assembly: Start=%g s, End=%g s, K_global_per_acq=%d samples.\n',min_global_start_time,max_global_end_time,K_global_per_acq);
H_assembled=zeros(K_global_per_acq*R_acquisitions,N_pixels); current_row_offset=0;
for r_acq=1:R_acquisitions
    hhp_current=all_hhp_data{r_acq}; start_time_current=all_start_times(r_acq); K_current=all_K_values(r_acq);
    if K_current==0||isempty(hhp_current); current_row_offset=current_row_offset+K_global_per_acq; continue; end
    t_current_acq_axis=start_time_current+(0:(K_current-1))/fs; hhp_aligned_r=zeros(K_global_per_acq,N_pixels);
    for px_col=1:N_pixels; if ~isempty(hhp_current)&&size(hhp_current,2)>=px_col; if length(t_current_acq_axis)>1&&issorted(t_current_acq_axis); hhp_aligned_r(:,px_col)=interp1(t_current_acq_axis,hhp_current(:,px_col),t_common_axis,'linear',0); elseif isscalar(t_current_acq_axis) && K_global_per_acq >=1 ; idx_match=find(abs(t_common_axis-t_current_acq_axis)<(0.5/fs),1); if ~isempty(idx_match); hhp_aligned_r(idx_match,px_col)=hhp_current(1,px_col); end; end; end; end
    H_assembled(current_row_offset+(1:K_global_per_acq),:)=hhp_aligned_r; current_row_offset=current_row_offset+K_global_per_acq;
end
H=H_assembled;
M=size(H,1);
N=N_pixels; 
fprintf('Final Assembled H Matrix: %d rows (M) x %d columns (N).\n',M,N); if M==0;error('H matrix zero rows.');end;fprintf('\n');

% --- Plot a few columns of H (Figure 3) ---
figure(3);clf;hold on;num_cols_to_plot_H=min(N,5);indices_to_plot_H=round(linspace(1,N,num_cols_to_plot_H));t_axis_H_plot_full_concat=(0:(M-1))/fs*1e6;
for n_idx_plot=1:length(indices_to_plot_H);col_idx=indices_to_plot_H(n_idx_plot);plot(t_axis_H_plot_full_concat,H(:,col_idx),'DisplayName',sprintf('H col Px %d',col_idx));end
hold off;xlabel('Overall Row Index (us)');ylabel('Amplitude');title('Figure 3: Columns of Assembled H');axis tight;grid on;legend('Location','best');set(gcf,'Color','w');
fprintf('\n');

% --- Generate Scene (v) ---
fprintf('\n--- Generating Scene (v) ---\n');
scene_xmin=min(x_coords_img);scene_xmax=max(x_coords_img);scene_zmin=min(z_coords_img);scene_zmax=max(z_coords_img);scene_z_range=scene_zmax-scene_zmin;
targets.x=[scene_xmin*0.3,0,scene_xmax*0.3];targets.y=zeros(size(targets.x));targets.z=[scene_zmin+0.25*scene_z_range,scene_zmin+0.5*scene_z_range,scene_zmin+0.75*scene_z_range];amplitudes=[1,0.9,1.1]';
scene_matrix=zeros(length(z_coords_img),length(x_coords_img));
for i=1:length(targets.x);[~,ix_scene]=min(abs(x_coords_img-targets.x(i)));[~,iz_scene]=min(abs(z_coords_img-targets.z(i)));if ix_scene>0&&ix_scene<=length(x_coords_img)&&iz_scene>0&&iz_scene<=length(z_coords_img);scene_matrix(iz_scene,ix_scene)=amplitudes(i);end;end
v_true_vector=scene_matrix(:);fprintf('Scene generated with %d targets.\n',length(targets.x));

figure(4);clf;subplot(1,2,1);imagesc(x_coords_img*1e3,z_coords_img*1e3,scene_matrix);axis image;colormap gray;colorbar;xlabel('x (mm)');ylabel('z (mm)');title('True Image v (matrix)');set(gca,'YDir','normal');
subplot(1,2,2);stem(v_true_vector);axis square tight;title('True Image v (vector)');sgtitle('Figure 4: True Image');set(gcf,'Color','w');set(gcf,'Position',[100,100,900,450]);fprintf('\n');

% --- Image Formation Model: u = H*v + n ---
fprintf('\n--- Simulating Measurement Acquisition (u = Hv + n) ---\n');
Hv_signal=H*v_true_vector;electronic_SNR=10^(electronic_SNR_db/10);signal_power_est=mean(Hv_signal(:).^2);
if isnan(signal_power_est)||signal_power_est<1e-30;warning('Signal power zero/NaN (%g). Small fallback noise.',signal_power_est);noise_sigma=1e-9*max(abs(H(:)));if abs(noise_sigma)<eps;noise_sigma=1e-15;end;actual_SNR_db=-Inf;
else;noise_variance=signal_power_est/electronic_SNR;noise_sigma=sqrt(noise_variance);actual_SNR_db=10*log10(signal_power_est/noise_variance);end
n_noise_vec=noise_sigma*randn(size(Hv_signal));u_measured_signal=Hv_signal+n_noise_vec;
fprintf('Noise: Target SNR %.1f dB. Actual SNR %.1f dB. Sigma: %g\n',electronic_SNR_db,actual_SNR_db,noise_sigma);

figure(6);clf;t_ax_meas=(0:(M-1))/fs*1e6;plot_lim_samp=min(M,round(fs*max_global_end_time*1.1));if plot_lim_samp==0&&M>0;plot_lim_samp=M;end
subplot(1,3,1);if ~isempty(Hv_signal)&&plot_lim_samp>0;plot(t_ax_meas(1:plot_lim_samp),Hv_signal(1:plot_lim_samp));end;axis tight;grid on;xlabel('Time (us)');title('Hv - Ideal');yl=ylim;
subplot(1,3,2);if ~isempty(n_noise_vec)&&plot_lim_samp>0;plot(t_ax_meas(1:plot_lim_samp),n_noise_vec(1:plot_lim_samp));end;axis tight;grid on;xlabel('Time (us)');title(sprintf('Noise (SNR~%.1fdB)',actual_SNR_db));if all(isfinite(yl));ylim(yl);end
subplot(1,3,3);if ~isempty(u_measured_signal)&&plot_lim_samp>0;plot(t_ax_meas(1:plot_lim_samp),u_measured_signal(1:plot_lim_samp));end;axis tight;grid on;xlabel('Time (us)');title('u (Meas+Noise)');if all(isfinite(yl));ylim(yl);end
sgtitle(sprintf('Figure 6: Measurement Signals (Compression=%.2f)',N/M));set(gcf,'Color','w');fprintf('\n');

%% Reconstruction Algorithms
fprintf('\n--- Starting Reconstruction ---\n');
A_matrix=H;
At_matrix=transpose(A_matrix);
b_vector=u_measured_signal;
I_true_matrix=scene_matrix;
v_true_vec_norm=v_true_vector;
if max(abs(v_true_vec_norm))>0
    v_true_vec_norm=v_true_vec_norm./max(abs(v_true_vec_norm));
end
imageResolution=size(I_true_matrix);
Afun=@(x_vec)A_matrix*x_vec;
Atfun=@(y_vec)At_matrix*y_vec;
use_pcg_regularization = false; 
alpha_pcg_reg = 1e-20; % PCG defaults

% --- Least Norm solution with PCG ---
fprintf('  Calculating Least Norm (PCG) Solution...\n'); tic;
AAtfun_pcg=@(y_vec)Afun(Atfun(y_vec));
fprintf('    PCG Parameters: maxIters=%d, tol=%g, Regularization: %s (alpha=%g)\n', maxItersCG_main, tolCG_main, mat2str(use_pcg_regularization), alpha_pcg_reg);
AAtfun_pcg_eff = AAtfun_pcg;
if use_pcg_regularization
    fprintf('    Using PCG with Tikhonov Regularization (alpha=%g)\n',alpha_pcg_reg);AAtfun_pcg_eff = @(y_vec)AAtfun_pcg(y_vec)+alpha_pcg_reg*y_vec;
end
x_pcg_img_vec=zeros(N,1);
iter_cg_main=0;
flag_cg_main=-1;
relres_cg_main=NaN;
if M==0||N==0||isempty(b_vector)||norm(b_vector(:))==0 
    fprintf('    Skipping PCG.\n');
else
    [y_sol,flag_cg_main,relres_cg_main,iter_cg_main]=pcg(AAtfun_pcg_eff,b_vector(:),tolCG_main,maxItersCG_main);
    fprintf('    PCG Finished: Flag=%d, Iterations=%d/%d, Relative Residual=%g\n',flag_cg_main,iter_cg_main,maxItersCG_main,relres_cg_main);
    if flag_cg_main~=0;warning('PCG did not converge as expected (Flag=%d).',flag_cg_main);
    end
    x_pcg_img_vec=Atfun(y_sol);
end
fprintf('    DEBUG PCG: Min raw x_pcg_img_vec: %g, Max raw: %g\n', min(x_pcg_img_vec(:)), max(x_pcg_img_vec(:)));
x_pcg_proc=real(x_pcg_img_vec);min_v=min(x_pcg_proc(:));max_v=max(x_pcg_proc(:));
if (max_v-min_v)>eps;x_pcg_norm=(x_pcg_proc-min_v)/(max_v-min_v);else;fprintf('    DEBUG PCG: Raw dynamic range too small. Normalized image set to zeros.\n');x_pcg_norm=zeros(N,1);end
x_pcg_norm(isnan(x_pcg_norm))=0;x_pcg2D=reshape(x_pcg_norm,imageResolution);
MSE_pcg=mean((x_pcg_norm-v_true_vec_norm).^2);PSNR_pcg=10*log10(1/MSE_pcg);runtime_pcg=toc;
fprintf('  PCG: MSE=%g, PSNR=%.2f dB, Iters Used=%d, Time=%.2f s\n',MSE_pcg,PSNR_pcg,iter_cg_main,runtime_pcg);fprintf('\n');

% --- Least Norm solution with Pinv ---
fprintf('  Calculating Least Norm (Pseudo-inverse) Solution...\n');tic; tolPinv=0;
x_pinv_img_vec=zeros(N,1);
pinv_method_desc='Skipped';
if M>0&&N>0&&~isempty(A_matrix)&&~isempty(b_vector); try if M<N; fprintf('    Pinv: M<N method.\n');AA_inv=pinv(A_matrix*At_matrix,tolPinv);x_pinv_img_vec=At_matrix*(AA_inv*b_vector(:));pinv_method_desc='At*pinv(A*At)*b'; else; fprintf('    Pinv: M>=N method.\n');A_pinv=pinv(A_matrix,tolPinv);x_pinv_img_vec=A_pinv*b_vector(:);pinv_method_desc='pinv(A)*b';end; catch ME; fprintf('    Pinv failed: %s\n',ME.message);pinv_method_desc='Failed';end;else;fprintf('    Skipping Pinv.\n');end
fprintf('    DEBUG Pinv: MinRaw=%g,MaxRaw=%g,NaNs=%d,Infs=%d\n',min(x_pinv_img_vec(:)),max(x_pinv_img_vec(:)),nnz(isnan(x_pinv_img_vec)),nnz(isinf(x_pinv_img_vec)));
x_pinv_proc=real(x_pinv_img_vec);min_v=min(x_pinv_proc(:));max_v=max(x_pinv_proc(:));
if(max_v-min_v)>eps;x_pinv_norm=(x_pinv_proc-min_v)/(max_v-min_v);else;x_pinv_norm=zeros(N,1);end;x_pinv_norm(isnan(x_pinv_norm))=0;
x_pinv2D=reshape(x_pinv_norm,imageResolution);MSE_pinv=mean((x_pinv_norm-v_true_vec_norm).^2);PSNR_pinv=10*log10(1/MSE_pinv);runtime_pinv=toc;
fprintf('  Pinv (%s): MSE=%g, PSNR=%.2f dB, Time=%.2f s\n',pinv_method_desc,MSE_pinv,PSNR_pinv,runtime_pinv);fprintf('\n');

% --- Plot Least Norm Results (Figure 7) ---
figure(7);clf;subplot(1,3,1);imagesc(x_coords_img*1e3,z_coords_img*1e3,I_true_matrix);axis image;colormap gray;colorbar;set(gca,'YDir','normal');title('Ground Truth');xlabel('x (mm)');ylabel('z (mm)');
subplot(1,3,2);imagesc(x_coords_img*1e3,z_coords_img*1e3,x_pcg2D);axis image;colormap gray;colorbar;set(gca,'YDir','normal');clim([0 1]);title(sprintf('Least Norm (PCG)\nMaxIter=%d, Tol=%g, UsedIter=%d\nRuntime=%.1fs, PSNR=%.2fdB',maxItersCG_main,tolCG_main,iter_cg_main,runtime_pcg,PSNR_pcg));xlabel('x (mm)');
subplot(1,3,3);imagesc(x_coords_img*1e3,z_coords_img*1e3,x_pinv2D);axis image;colormap gray;colorbar;set(gca,'YDir','normal');clim([0 1]);title(sprintf('Least Norm (Pinv)\nTol=%g, Runtime=%.1fs\nPSNR=%.2fdB',tolPinv,runtime_pinv,PSNR_pinv));xlabel('x (mm)');
sgtitle(sprintf('Figure 7: Least Norm Reconstructions (R_{acq}=%d, SNR~%.1fdB)',R_acquisitions,actual_SNR_db));set(gcf,'Color','w');set(gcf,'Position',[150,150,1300,450]);

% --- ADMM with TV regularization ---
fprintf('  Calculating ADMM-TV Solution...\n');
fprintf('    ADMM Parameters: numIters=%d, rho=%g, lambda_TV=%g\n', numItersADMM, rho_admm, lambda_tv_reg);
H_norm_factor=max(abs(A_matrix(:))); if H_norm_factor<eps;H_norm_factor=1;end; A_admm=A_matrix./H_norm_factor; At_admm=transpose(A_admm);b_admm_vec=b_vector(:)/H_norm_factor;noise_sigma_admm=noise_sigma/H_norm_factor;Afun_admm=@(x)A_admm*x(:); Atfun_admm_img=@(y)reshape(At_admm*y,imageResolution); AtAfun_admm_img=@(x)Atfun_admm_img(Afun_admm(x));[Dx_sparse,Dy_sparse]=createDifferenceOperators(imageResolution); opDx_tv=@(x)[Dx_sparse*x(:),Dy_sparse*x(:)];opDtx_tv=@(v)reshape(Dx_sparse'*v(:,1)+Dy_sparse'*v(:,2),imageResolution);opDtDx_tv=@(x)opDtx_tv(opDx_tv(x));
x_admm_img_iter=zeros(imageResolution);z_admm_grad_iter=zeros([prod(imageResolution) 2]);u_admm_dual_iter=zeros([prod(imageResolution) 2]);Hfun_pcg_admm=@(x_vec)reshape(AtAfun_admm_img(reshape(x_vec,imageResolution))+rho_admm.*opDtDx_tv(reshape(x_vec,imageResolution)),[prod(imageResolution) 1]);PSNR_admm_iters=zeros([numItersADMM 1]);residuals_admm_iters=zeros([numItersADMM 1]);
figure(8);clf;set(gcf,'Position',[200,200,1000,700],'Color','w');tic;
for k_admm=1:numItersADMM; figure(8); % Ensure Figure 8 is current
    v_upd=z_admm_grad_iter-u_admm_dual_iter;bb_upd=Atfun_admm_img(b_admm_vec)+rho_admm*opDtx_tv(v_upd);[x_vec_new,~,~,~]=pcg(Hfun_pcg_admm,bb_upd(:),noise_sigma_admm,25,[],[],x_admm_img_iter(:));x_admm_img_iter=reshape(x_vec_new,imageResolution);kap=lambda_tv_reg/rho_admm;v_z_upd=opDx_tv(x_admm_img_iter)+u_admm_dual_iter;v_norm=sqrt(sum(v_z_upd.^2,2));v_norm(v_norm<eps)=eps;shr=max(0,1-kap./v_norm);z_admm_grad_iter=v_z_upd.*shr;u_admm_dual_iter=u_admm_dual_iter+opDx_tv(x_admm_img_iter)-z_admm_grad_iter;x_scl=real(x_admm_img_iter(:));min_x=min(x_scl);max_x=max(x_scl);if(max_x-min_x)>eps;x_scl=(x_scl-min_x)/(max_x-min_x);else;x_scl=zeros(size(x_scl));end;MSE_curr=mean((x_scl-v_true_vec_norm).^2);PSNR_admm_iters(k_admm)=10*log10(1/MSE_curr);r1=b_admm_vec-Afun_admm(x_admm_img_iter);r2=opDx_tv(x_admm_img_iter);tv_n=sum(sqrt(sum(r2.^2,2)));residuals_admm_iters(k_admm)=0.5*sum(r1(:).^2)+lambda_tv_reg*tv_n;
    subplot(2,3,[1 4]);imagesc(x_coords_img*1e3,z_coords_img*1e3,I_true_matrix);axis image;colormap(gca,gray);colorbar;set(gca,'YDir','normal');title('Target');xlabel('x(mm)');ylabel('z(mm)');subplot(2,3,[2 5]);imagesc(x_coords_img*1e3,z_coords_img*1e3,reshape(x_scl,imageResolution));axis image;colormap(gca,gray);colorbar;set(gca,'YDir','normal');title(sprintf('\\lambda=%.1e,\\rho=%.1f\nPSNR=%.2fdB,It %d',lambda_tv_reg,rho_admm,PSNR_admm_iters(k_admm),k_admm));xlabel('x(mm)');subplot(2,3,3);plot(1:k_admm,PSNR_admm_iters(1:k_admm),'m-','LineWidth',2);title('PSNR/Iter');xlabel('Iter');ylabel('PSNR(dB)');grid on;axis tight;if k_admm>1;yl=ylim;if diff(yl)<1;yl(2)=yl(1)+1;end;ylim(yl);end;subplot(2,3,6);plot(1:k_admm,log10(residuals_admm_iters(1:k_admm)),'-','LineWidth',2);title('log10(Obj)/Iter');xlabel('Iter');ylabel('log10(Val)');grid on;axis tight;if k_admm>1;yl=ylim;if diff(yl)<0.1;yl(2)=yl(1)+0.1;end;ylim(yl);end;drawnow;
end
runtime_ADMM_total=toc;fprintf('  ADMM-TV: Final PSNR=%.2f dB, Time=%.2fs\n',PSNR_admm_iters(end),runtime_ADMM_total);sgtitle(sprintf('Figure 8: ADMM TV (R_{acq}=%d, SNR~%.1fdB)',R_acquisitions,actual_SNR_db));fprintf('\n');

%% End Field II Simulation
if exist('field_end','file')==2;field_end;disp('Field II ended.');else;disp('field_end function not found.');end;fprintf('\n');

%% Helper function
function [Dx,Dy]=createDifferenceOperators(imageSize)
    rows=imageSize(1);cols=imageSize(2);N_img_pixels=rows*cols;
    Dx=spdiags([-ones(N_img_pixels,1),ones(N_img_pixels,1)],[0,rows],N_img_pixels,N_img_pixels);
    last_col_indices_mask=false(N_img_pixels,1);last_col_indices_mask((cols-1)*rows+1:cols*rows)=true;Dx(last_col_indices_mask,:)=0;
    Dy=spdiags([-ones(N_img_pixels,1),ones(N_img_pixels,1)],[0,1],N_img_pixels,N_img_pixels);
    last_row_indices_mask=false(N_img_pixels,1);last_row_indices_mask(rows:rows:N_img_pixels)=true;Dy(last_row_indices_mask,:)=0;
end