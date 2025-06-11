% Advanced_pMUT_CS_Recon_v5.3_Final.m
%
% Description:
% This version contains the final fix. It keeps the high-energy burst
% excitation from v5.2 and corrects the logic for noise calculation.
%
% --- CRITICAL CHANGE ---
% 1. The threshold for detecting signal power is lowered from 1e-20 to 1e-30,
%    allowing the script to correctly scale the noise to our ~10^-12 signal.
%
clearvars; clc; close all;
% --- Initialize Field II ---
field_init(-1);
% --- Transducer Parameters (Temporal) ---
fc_nominal = 5.8e4; % Nominal Center Frequency [Hz] (Used for lambda reference)
c = 1540;           % Speed of Sound [m/s]
lambda = c/fc_nominal; % Wavelength [m] (based on nominal fc)
fs = 2e6;           % Simulation Sampling Frequency [Hz]. High fixed value for chirp.
fprintf('Simulation sampling frequency fs = %g MHz\n', fs/1e6);
% --- Transducer Parameters (Spatial - for 3-pMUT Array) ---
pMUT_width = 1.5*lambda;    % Width of each pMUT element
pMUT_height = pMUT_width;   % Height of each pMUT element
kerf = 0.5*lambda;          % Kerf between pMUT elements (if they were in a dense grid)
d_spacing = 5*lambda;       % Characteristic spacing for the triangular pMUT array
% Define desired positions for a 3-element triangular array
pos1 = [d_spacing,0,0];
pos2 = [d_spacing*cos(2*pi/3), d_spacing*sin(2*pi/3),0];
pos3 = [d_spacing*cos(4*pi/3), d_spacing*sin(4*pi/3),0];
desired_positions = [pos1; pos2; pos3];
num_active_intended = size(desired_positions, 1);
fprintf('Intending to use %d pMUTs in a triangular configuration.\n', num_active_intended);
% Define a Field II grid large enough to encompass these desired positions
array_span_x = 2.5*d_spacing; % Span of the encompassing grid
array_span_y = 2.5*d_spacing;
element_width_grid = pMUT_width;  % Grid element size matches pMUT size
element_height_grid = pMUT_height;
num_x_grid = ceil(array_span_x / (element_width_grid + kerf));
num_y_grid = ceil(array_span_y / (element_height_grid + kerf));
if mod(num_x_grid, 2) == 0; num_x_grid = num_x_grid + 1; end % Ensure odd for centering if needed
if mod(num_y_grid, 2) == 0; num_y_grid = num_y_grid + 1; end
% --- Voxels (Imaging Grid) ---
grid_width = 2*lambda;         % Lateral width of the imaging region
grid_depth_start = 12*lambda;   % Starting depth of the imaging region
grid_depth_end = 40*lambda;     % Ending depth of the imaging region
grid_step = lambda/2;           % Pixel size (resolution)
x_coords_img = -grid_width/2 : grid_step : grid_width/2;
z_coords_img = grid_depth_start : grid_step : grid_depth_end;
[X_mesh, Z_mesh] = meshgrid(x_coords_img, z_coords_img);
Y_mesh = zeros(size(X_mesh)); % For 2D imaging in x-z plane
N_pixels = numel(X_mesh); % Total number of pixels in the image
fprintf('Wavelength lambda (nominal) = %g mm\n', lambda*1e3);
fprintf('Imaging Grid: %d pixels (%d axial x %d lateral).\n', N_pixels, length(z_coords_img), length(x_coords_img));
% --- Set Field Parameters ---
set_field('fs', fs);
set_field('c', c);
% --- Define Impulse Response and Global Excitation ---
impulseResponse = 1;
fprintf('Set impulse response to a delta function (broadband transducer model).\n');
excitationPulse = 1; % Global excitation (largely overridden by ele_waveform)
% --- Define pMUT Aperture (Tx and Rx are the same) ---
physical_element_centers = zeros(num_x_grid * num_y_grid, 3);
element_no_grid_map = 0;
center_offset_x = (num_x_grid - 1)/2 * (element_width_grid + kerf);
center_offset_y = (num_y_grid - 1)/2 * (element_height_grid + kerf);
for iy = 1:num_y_grid
    y_pos_el = (iy-1)*(element_height_grid + kerf) - center_offset_y;
    for ix = 1:num_x_grid
        x_pos_el = (ix-1)*(element_width_grid + kerf) - center_offset_x;
        element_no_grid_map = element_no_grid_map + 1;
        physical_element_centers(element_no_grid_map, :) = [x_pos_el, y_pos_el, 0];
    end
end
active_indices_linear = zeros(num_active_intended, 1);
for i = 1:num_active_intended
    distances = sqrt(sum((physical_element_centers - desired_positions(i,:)).^2, 2));
    [~, min_idx] = min(distances);
    active_indices_linear(i) = min_idx;
end
active_indices_linear = unique(active_indices_linear);
num_active = length(active_indices_linear);
if num_active ~= num_active_intended; warning('Could not map to %d unique elements (found %d).', num_active_intended, num_active); if num_active == 0; error('No active elements found!'); end; end
fprintf('Successfully mapped %d unique active pMUTs.\n', num_active);
enabled_matrix = zeros(num_y_grid, num_x_grid);
[row_indices, col_indices] = ind2sub([num_y_grid, num_x_grid], active_indices_linear);
for i = 1:num_active; enabled_matrix(row_indices(i), col_indices(i)) = 1; end
nominal_focus_point = [0 0 100];
pMUT_Aperture = xdc_2d_array(num_x_grid, num_y_grid, element_width_grid, element_height_grid, kerf, kerf, enabled_matrix, 1, 1, nominal_focus_point);
xdc_impulse(pMUT_Aperture, impulseResponse);
xdc_excitation(pMUT_Aperture, excitationPulse);
% --- Plot pMUT Array Geometry (Figure 1) ---
figure(1); clf;
[rows_active_plot, cols_active_plot] = find(enabled_matrix);
active_centers_plot = zeros(num_active, 3);
for i = 1:num_active; iy_plot=rows_active_plot(i);ix_plot=cols_active_plot(i);y_pos_plot=(iy_plot-1)*(element_height_grid+kerf)-center_offset_y;x_pos_plot=(ix_plot-1)*(element_width_grid+kerf)-center_offset_x;active_centers_plot(i,:)=[x_pos_plot,y_pos_plot,0];end
plot3(physical_element_centers(:,1)*1e3,physical_element_centers(:,2)*1e3,physical_element_centers(:,3)*1e3,'k.','MarkerSize',1,'DisplayName','Grid Elements'); hold on;
plot3(active_centers_plot(:,1)*1e3,active_centers_plot(:,2)*1e3,active_centers_plot(:,3)*1e3,'ro','MarkerSize',8,'MarkerFaceColor','r','DisplayName','Active pMUTs');
xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)');title(sprintf('pMUT Array Geometry (%d Active Elements)',num_active));
axis equal; grid on; view(2); legend('Location','best'); set(gcf,'Color','w');

%% Generate Coded Acquisitions and Measure H components
R_acquisitions = 20;
fprintf('\n--- Generating H Matrix Components from %d Coded Acquisitions ---\n', R_acquisitions);
hydrophone_positions_img = [X_mesh(:), Y_mesh(:), Z_mesh(:)];

% --- High-Energy Burst Excitation ---
f_start_chirp = 10e3;    % 10 kHz
f_end_chirp = 200e3;   % 200 kHz
burst_duration = 0.1e-3; % 0.1 ms (100 us) burst duration
excitation_amplitude = 500; % KEY PARAMETER: Numerical amplitude boost
fprintf('Synthetic Burst: %g-%g kHz, duration %g us, amplitude %g\n', f_start_chirp/1e3, f_end_chirp/1e3, burst_duration*1e3, excitation_amplitude);
t_chirp_vec = 0:1/fs:burst_duration;
synth_burst_base = chirp(t_chirp_vec, f_start_chirp, t_chirp_vec(end), f_end_chirp, 'linear');
synth_burst_windowed = synth_burst_base .* tukeywin(length(t_chirp_vec), 0.25)';
base_chirp_to_modulate = synth_burst_windowed(:)' * excitation_amplitude;

all_hhp_data=cell(R_acquisitions,1); all_start_times=zeros(R_acquisitions,1); all_K_values=zeros(R_acquisitions,1);
rng('default');

% --- Plot Example Transmitted Field (Figure 2 adaptation) ---
figure(2); clf;
if num_active > 0
    ele_waveform(pMUT_Aperture, 1, base_chirp_to_modulate);
    for p_idx_other=2:num_active; ele_waveform(pMUT_Aperture,p_idx_other,0); end
    first_pmut_pos=active_centers_plot(1,:); field_calc_x=first_pmut_pos(1); field_calc_y=first_pmut_pos(2);
    field_calc_z=(grid_depth_start:lambda/4:grid_depth_end);
    field_calc_points=[repmat(field_calc_x,length(field_calc_z),1),repmat(field_calc_y,length(field_calc_z),1),field_calc_z(:)];
    if ~isempty(field_calc_points)
        [hp_example,~]=calc_hp(pMUT_Aperture,field_calc_points);
        hp_env_example=abs(hilbert(hp_example));
        imagesc(t_chirp_vec*1e6,field_calc_z*1e3,hp_env_example');
        xlabel('Time (us)');ylabel('Depth (mm)');title('Example Transmitted Field (Envelope from 1st pMUT)');colorbar;
    else; title('Example Transmitted Field (Skipped - No field points)'); 
    end
else; title('Example Transmitted Field (Skipped - No active pMUTs)'); 
end
sgtitle('Figure 2: Example Transmitted Field'); set(gcf,'Color','w');

% --- Loop for each coded acquisition ---
for p_idx=1:num_active; ele_waveform(pMUT_Aperture,p_idx,0); end % Reset
for r_acq=1:R_acquisitions
    fprintf('Calculating Acquisition %d/%d...\n',r_acq,R_acquisitions);
    for p_idx=1:num_active
        current_rand_phase=2*pi*rand();
        analytic_signal=hilbert(base_chirp_to_modulate);
        phase_shifted_analytic=analytic_signal*exp(1i*current_rand_phase);
        current_pmut_waveform=real(phase_shifted_analytic);
        ele_waveform(pMUT_Aperture,p_idx,current_pmut_waveform);
    end
    [hhp_r,start_time_r]=calc_hhp(pMUT_Aperture,pMUT_Aperture,hydrophone_positions_img);
    all_hhp_data{r_acq}=hhp_r; all_start_times(r_acq)=start_time_r; all_K_values(r_acq)=size(hhp_r,1);
    if ~isempty(hhp_r); fprintf('    Max abs value in hhp_r: %g (K_r=%d)\n',max(abs(hhp_r(:))),size(hhp_r,1));
    else; fprintf('    WARNING: hhp_r is empty for acquisition %d!\n',r_acq); end
end

%% Assemble Full H Matrix with Interpolation and Alignment
fprintf('\n--- Assembling Full H Matrix ---\n');
all_end_times=zeros(R_acquisitions,1); for r_idx=1:R_acquisitions; if all_K_values(r_idx)>0; all_end_times(r_idx)=all_start_times(r_idx)+(all_K_values(r_idx)-1)/fs; end; end
min_global_start_time=min(all_start_times); max_global_end_time=max(all_end_times);
if min_global_start_time>=max_global_end_time; if any(all_K_values>0); max_global_end_time=min_global_start_time+(max(all_K_values(all_K_values>0))-1)/fs; else; min_global_start_time=0;max_global_end_time=100/fs; warning('No data in any acquisition, H will be all zeros.'); end; end
t_common_axis=min_global_start_time:1/fs:max_global_end_time; K_global_per_acq=length(t_common_axis);
fprintf('Global Time Window for H assembly: Start=%g s, End=%g s, K_global_per_acq=%d samples.\n',min_global_start_time,max_global_end_time,K_global_per_acq);
H_assembled=zeros(K_global_per_acq*R_acquisitions,N_pixels); current_row_offset=0;
for r_acq=1:R_acquisitions
    hhp_current=all_hhp_data{r_acq}; start_time_current=all_start_times(r_acq); K_current=all_K_values(r_acq);
    if K_current==0||isempty(hhp_current); current_row_offset=current_row_offset+K_global_per_acq; continue; end
    t_current_acq_axis=start_time_current+(0:(K_current-1))/fs;
    hhp_aligned_r=zeros(K_global_per_acq,N_pixels);
    for px_col=1:N_pixels
        if ~isempty(hhp_current)&&size(hhp_current,2)>=px_col; if length(t_current_acq_axis)>1&&issorted(t_current_acq_axis); hhp_aligned_r(:,px_col)=interp1(t_current_acq_axis,hhp_current(:,px_col),t_common_axis,'linear',0); elseif isscalar(t_current_acq_axis); idx_match=find(abs(t_common_axis-t_current_acq_axis)<(0.5/fs),1); if ~isempty(idx_match); hhp_aligned_r(idx_match,px_col)=hhp_current(1,px_col); end; end; end
    end
    H_assembled(current_row_offset+(1:K_global_per_acq),:)=hhp_aligned_r; current_row_offset=current_row_offset+K_global_per_acq;
end
H=H_assembled; M=size(H,1); N=N_pixels;
fprintf('Final Assembled H Matrix: %d rows (M) x %d columns (N).\n', M, N);
if M==0; error('H matrix has zero rows. Cannot proceed.'); end

% --- Plot a few columns of H (Figure 3 adaptation) ---
figure(3); clf; hold on; num_cols_to_plot_H=min(N,5);
indices_to_plot_H=round(linspace(1,N,num_cols_to_plot_H));
t_axis_H_plot_full_concat=(0:(M-1))/fs*1e6;
for n_idx_plot=1:length(indices_to_plot_H); col_idx=indices_to_plot_H(n_idx_plot); plot(t_axis_H_plot_full_concat,H(:,col_idx),'DisplayName',sprintf('H column for Pixel %d',col_idx)); end
hold off; xlabel('Overall Row Index (Effective Concatenated Time - us)'); ylabel('Amplitude'); title('Figure 3: Columns of Assembled H Matrix');
axis tight;grid on;legend('Location','best');set(gcf,'Color','w');

%% Generate a Scene to be Imaged (v)
fprintf('\n--- Generating Scene (v) ---\n');
scene_xmin=min(x_coords_img);scene_xmax=max(x_coords_img);scene_zmin=min(z_coords_img);scene_zmax=max(z_coords_img);scene_z_range=scene_zmax-scene_zmin;
targets.x=[scene_xmin*0.5,0,scene_xmax*0.5,scene_xmin*0.2];targets.y=zeros(size(targets.x));
targets.z=[scene_zmin+0.2*scene_z_range,scene_zmin+0.5*scene_z_range,scene_zmin+0.8*scene_z_range,scene_zmin+0.4*scene_z_range];
amplitudes=[1,0.9,1.1,0.7]';
scene_matrix=zeros(length(z_coords_img),length(x_coords_img));
for i=1:length(targets.x);[~,ix_scene]=min(abs(x_coords_img-targets.x(i)));[~,iz_scene]=min(abs(z_coords_img-targets.z(i)));if ix_scene>0&&ix_scene<=length(x_coords_img)&&iz_scene>0&&iz_scene<=length(z_coords_img);scene_matrix(iz_scene,ix_scene)=amplitudes(i);end;end
v_true_vector=scene_matrix(:);
fprintf('Scene generated with %d targets.\n',length(targets.x));
figure(4);clf;subplot(1,2,1);imagesc(x_coords_img*1e3,z_coords_img*1e3,scene_matrix);axis image;colormap gray;colorbar;xlabel('x (mm)');ylabel('z (mm)');title('True Image v (matrix form)');set(gca,'YDir','normal');
subplot(1,2,2);stem(v_true_vector);axis square tight;title('True Image v (vector form)');xlabel('Pixel Index');ylabel('Amplitude');sgtitle('Figure 4: True Image (Ground Truth)');set(gcf,'Color','w');

%% Image Formation Model: u = H*v + n
fprintf('\n--- Simulating Measurement Acquisition (u = Hv + n) ---\n');
Hv_signal=H*v_true_vector;
rng('shuffle');electronic_SNR_db=90;electronic_SNR=10^(electronic_SNR_db/10);
signal_power_est=mean(Hv_signal(:).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE FINAL FIX: Lowering the signal power threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(signal_power_est) || signal_power_est < 1e-30 % OLD THRESHOLD: 1e-20
    warning('Signal power is effectively zero or NaN (%g). Using small default noise sigma.', signal_power_est);
    noise_sigma = 1e-9; % Fallback noise sigma
    actual_SNR_db = -Inf;
else
    % This block should now execute correctly
    noise_variance = signal_power_est / electronic_SNR;
    noise_sigma = sqrt(noise_variance);
    actual_SNR_db = 10*log10(signal_power_est / noise_variance);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_noise_vec=noise_sigma*randn(size(Hv_signal));u_measured_signal=Hv_signal+n_noise_vec;
fprintf('Added Gaussian noise. Target SNR: %.1f dB. Actual SNR: %.1f dB. Noise sigma: %g\n',electronic_SNR_db,actual_SNR_db,noise_sigma);
figure(6);clf;t_axis_measurements_plot=(0:(M-1))/fs*1e6;
plot_limit_samples=min(M,round(fs*max_global_end_time*1.1));if plot_limit_samples==0&&M>0;plot_limit_samples=M;end
subplot(1,3,1);if ~isempty(Hv_signal)&&plot_limit_samples>0;plot(t_axis_measurements_plot(1:plot_limit_samples),Hv_signal(1:plot_limit_samples));else;plot([]);end;axis tight;grid on;xlabel('Time (µs)');title('Measurements (Hv - Ideal)');ylabel('Amplitude');yl_hv_plot=ylim;
subplot(1,3,2);if ~isempty(n_noise_vec)&&plot_limit_samples>0;plot(t_axis_measurements_plot(1:plot_limit_samples),n_noise_vec(1:plot_limit_samples));else;plot([]);end;axis tight;grid on;xlabel('Time (µs)');title(sprintf('Noise (SNR \\approx %.1f dB)',actual_SNR_db));ylabel('Amplitude');if ~isempty(yl_hv_plot)&&all(isfinite(yl_hv_plot));ylim(yl_hv_plot);end
subplot(1,3,3);if ~isempty(u_measured_signal)&&plot_limit_samples>0;plot(t_axis_measurements_plot(1:plot_limit_samples),u_measured_signal(1:plot_limit_samples));else;plot([]);end;axis tight;grid on;xlabel('Time (µs)');title('Measurements + Noise (u)');ylabel('Amplitude');if ~isempty(yl_hv_plot)&&all(isfinite(yl_hv_plot));ylim(yl_hv_plot);end
sgtitle(sprintf('Figure 6: Measurement Signals (Compression=%.2f)',N/M));set(gcf,'Color','w');

%% Reconstruction Algorithms
fprintf('\n--- Starting Reconstruction ---\n');
A_matrix=H;At_matrix=transpose(A_matrix);b_vector=u_measured_signal;I_true_matrix=scene_matrix;
v_true_vec_norm=v_true_vector;if max(abs(v_true_vec_norm))>0;v_true_vec_norm=v_true_vec_norm./max(abs(v_true_vec_norm));end
imageResolution=size(I_true_matrix);Afun=@(x_vec)A_matrix*x_vec;Atfun=@(y_vec)At_matrix*y_vec;
fprintf('  Calculating Least Norm (PCG) Solution...\n');tic;
AAtfun_pcg=@(y_intermediate_cg_vec)Afun(Atfun(y_intermediate_cg_vec));maxItersCG_main=200;tolCG_main=noise_sigma;
if M==0||N==0||isempty(b_vector)||norm(b_vector)==0;x_pcg_img_vec=zeros(N,1);iter_cg_main=0;flag_cg_main=-1;
else;[y_cg_solution_M_by_1,flag_cg_main,~,iter_cg_main]=pcg(AAtfun_pcg,b_vector(:),tolCG_main,maxItersCG_main);x_pcg_img_vec=Atfun(y_cg_solution_M_by_1);end
x_pcg_img_vec_processed=real(x_pcg_img_vec);min_val_pcg=min(x_pcg_img_vec_processed(:));max_val_pcg=max(x_pcg_img_vec_processed(:));
if (max_val_pcg-min_val_pcg)>eps;x_pcg_norm_vec=(x_pcg_img_vec_processed-min_val_pcg)/(max_val_pcg-min_val_pcg);else;x_pcg_norm_vec=zeros(N,1);end;x_pcg_norm_vec(isnan(x_pcg_norm_vec))=0;
x_pcg2D=reshape(x_pcg_norm_vec,imageResolution);MSE_pcg=mean((x_pcg_norm_vec-v_true_vec_norm).^2);PSNR_pcg=10*log10(1/MSE_pcg);runtime_pcg=toc;
fprintf('  PCG: MSE=%g, PSNR=%.2f dB, Iters=%d, Time=%.2f s\n',MSE_pcg,PSNR_pcg,iter_cg_main,runtime_pcg);
fprintf('  Calculating Least Norm (Pseudo-inverse) Solution...\n');tic;
tolPinv=1e-6;x_pinv_img_vec=zeros(N,1);pinv_method_desc='Skipped';
if M>0&&N>0&&~isempty(A_matrix)&&~isempty(b_vector);try if M<N;AAt_pinv_calc=A_matrix*At_matrix;AAt_inv=pinv(AAt_pinv_calc,tolPinv);x_pinv_img_vec=At_matrix*(AAt_inv*b_vector(:));pinv_method_desc='At*pinv(A*At)*b';else;A_pinv_calc=pinv(A_matrix,tolPinv);x_pinv_img_vec=A_pinv_calc*b_vector(:);pinv_method_desc='pinv(A)*b';end;catch ME;fprintf('    Pinv failed: %s\n',ME.message);pinv_method_desc='Failed';end;end
x_pinv_img_vec_processed=real(x_pinv_img_vec);min_val_pinv=min(x_pinv_img_vec_processed(:));max_val_pinv=max(x_pinv_img_vec_processed(:));
if (max_val_pinv-min_val_pinv)>eps;x_pinv_norm_vec=(x_pinv_img_vec_processed-min_val_pinv)/(max_val_pinv-min_val_pinv);else;x_pinv_norm_vec=zeros(N,1);end;x_pinv_norm_vec(isnan(x_pinv_norm_vec))=0;
x_pinv2D=reshape(x_pinv_norm_vec,imageResolution);MSE_pinv=mean((x_pinv_norm_vec-v_true_vec_norm).^2);PSNR_pinv=10*log10(1/MSE_pinv);runtime_pinv=toc;
fprintf('  Pinv (%s): MSE=%g, PSNR=%.2f dB, Time=%.2f s\n',pinv_method_desc,MSE_pinv,PSNR_pinv,runtime_pinv);
figure(7);clf;subplot(1,3,1);imagesc(x_coords_img*1e3,z_coords_img*1e3,I_true_matrix);axis image;colormap gray;colorbar;set(gca,'YDir','normal');title('Ground Truth');xlabel('x (mm)');ylabel('z (mm)');
subplot(1,3,2);imagesc(x_coords_img*1e3,z_coords_img*1e3,x_pcg2D);axis image;colormap gray;colorbar;set(gca,'YDir','normal');title(sprintf('Least Norm (PCG)\nPSNR=%.2fdB',PSNR_pcg));
subplot(1,3,3);imagesc(x_coords_img*1e3,z_coords_img*1e3,x_pinv2D);axis image;colormap gray;colorbar;set(gca,'YDir','normal');title(sprintf('Least Norm (Pinv)\nPSNR=%.2fdB',PSNR_pinv));
sgtitle('Figure 7: Least Norm Reconstructions');set(gcf,'Color','w');

fprintf('  Calculating ADMM-TV Solution...\n');
H_norm_factor=max(abs(A_matrix(:)));if H_norm_factor<eps;H_norm_factor=1;end;A_admm=A_matrix./H_norm_factor;At_admm=transpose(A_admm);
b_admm_vec=b_vector(:)/H_norm_factor;noise_sigma_admm=noise_sigma/H_norm_factor;
Afun_admm=@(x_img_vec)A_admm*x_img_vec;Atfun_admm_img=@(y_vec)reshape(At_admm*y_vec,imageResolution);
AtAfun_admm_img=@(x_img)Atfun_admm_img(Afun_admm(x_img(:)));[Dx_sparse,Dy_sparse]=createDifferenceOperators(imageResolution);
opDx_tv=@(x_img_tv)[Dx_sparse*x_img_tv(:),Dy_sparse*x_img_tv(:)];opDtx_tv=@(v_grad_vecs_tv)reshape(Dx_sparse'*v_grad_vecs_tv(:,1)+Dy_sparse'*v_grad_vecs_tv(:,2),imageResolution);
opDtDx_tv=@(x_img_tv)opDtx_tv(opDx_tv(x_img_tv));numItersADMM=25;rho_admm=10;lambda_tv_reg=0.01;
x_admm_img_iter=zeros(imageResolution);z_admm_grad_iter=zeros([prod(imageResolution) 2]);u_admm_dual_iter=zeros([prod(imageResolution) 2]);
Hfun_pcg_admm=@(x_img_vec_pcg)reshape(AtAfun_admm_img(reshape(x_img_vec_pcg,imageResolution))+rho_admm.*opDtDx_tv(reshape(x_img_vec_pcg,imageResolution)),[prod(imageResolution) 1]);
PSNR_admm_iters=zeros([numItersADMM 1]);residuals_admm_iters=zeros([numItersADMM 1]);
figure(8);clf;tic;
for k_admm=1:numItersADMM
    v_for_x_update=z_admm_grad_iter-u_admm_dual_iter;bb_for_x_update=Atfun_admm_img(b_admm_vec)+rho_admm*opDtx_tv(v_for_x_update);
    [x_admm_img_vec_new,~,~,~]=pcg(Hfun_pcg_admm,bb_for_x_update(:),noise_sigma_admm,25,[],[],x_admm_img_iter(:));
    x_admm_img_iter=reshape(x_admm_img_vec_new,imageResolution);
    kappa_tv_shrink=lambda_tv_reg/rho_admm;v_for_z_update=opDx_tv(x_admm_img_iter)+u_admm_dual_iter;
    v_grad_norm_z=sqrt(sum(v_for_z_update.^2,2));v_grad_norm_z(v_grad_norm_z<eps)=eps;
    shrinkage_factor_z=max(0,1-kappa_tv_shrink./v_grad_norm_z);z_admm_grad_iter=v_for_z_update.*shrinkage_factor_z;
    u_admm_dual_iter=u_admm_dual_iter+opDx_tv(x_admm_img_iter)-z_admm_grad_iter;
    x_admm_scaled_vec=real(x_admm_img_iter(:));min_x_admm_iter=min(x_admm_scaled_vec);max_x_admm_iter=max(x_admm_scaled_vec);
    if(max_x_admm_iter-min_x_admm_iter)>eps;x_admm_scaled_vec=(x_admm_scaled_vec-min_x_admm_iter)/(max_x_admm_iter-min_x_admm_iter);else;x_admm_scaled_vec=zeros(size(x_admm_scaled_vec));end
    MSE_current_admm=mean((x_admm_scaled_vec-v_true_vec_norm).^2);PSNR_admm_iters(k_admm)=10*log10(1/MSE_current_admm);
    r1_fidelity=b_admm_vec-Afun_admm(x_admm_img_iter(:));r2_tv_opDx_val=opDx_tv(x_admm_img_iter);
    tv_norm_val=sum(sqrt(sum(r2_tv_opDx_val.^2,2)));residuals_admm_iters(k_admm)=0.5*sum(r1_fidelity(:).^2)+lambda_tv_reg*tv_norm_val;
    figure(8)
    subplot(2,3,[1 4]);imagesc(x_coords_img*1e3,z_coords_img*1e3,I_true_matrix);axis image;colormap(gca,gray);colorbar;set(gca,'YDir','normal');title('Target Image');xlabel('x (mm)');ylabel('z (mm)');
    subplot(2,3,[2 5]);imagesc(x_coords_img*1e3,z_coords_img*1e3,reshape(x_admm_scaled_vec,imageResolution));axis image;colormap(gca,gray);colorbar;set(gca,'YDir','normal');title(sprintf('\\lambda=%.2e, \\rho=%.1f\nPSNR=%.2fdB, Iter %d',lambda_tv_reg,rho_admm,PSNR_admm_iters(k_admm),k_admm));xlabel('x (mm)');ylabel('z (mm)');
    subplot(2,3,3);plot(1:k_admm,PSNR_admm_iters(1:k_admm),'LineWidth',2,'color',[1 0 1]);title('PSNR vs Iteration');xlabel('Iteration');ylabel('PSNR (dB)');grid on;axis tight;
    subplot(2,3,6);plot(1:k_admm,log10(residuals_admm_iters(1:k_admm)),'LineWidth',2);title('log10(Objective) vs Iteration');xlabel('Iteration');ylabel('log10(Value)');grid on;axis tight;
    drawnow;
end
runtime_ADMM_total=toc;fprintf('  ADMM-TV complete: Final PSNR=%.2f dB, Total Time=%.2fs\n',PSNR_admm_iters(end),runtime_ADMM_total);
sgtitle(sprintf('Figure 8: ADMM TV (R_{acq}=%d, SNR~%.1fdB)',R_acquisitions,actual_SNR_db));set(gcf,'Color','w');

%% End Field II Simulation
if exist('field_end','file')==2;field_end;disp('Field II ended.');else;disp('field_end function not found.');end

%% Helper function
function [Dx,Dy]=createDifferenceOperators(imageSize)
rows=imageSize(1);cols=imageSize(2);N_img_pixels=rows*cols;
Dx=spdiags([-ones(N_img_pixels,1),ones(N_img_pixels,1)],[0,rows],N_img_pixels,N_img_pixels);
last_col_indices_mask=false(N_img_pixels,1);last_col_indices_mask((cols-1)*rows+1:cols*rows)=true;Dx(last_col_indices_mask,:)=0;
Dy=spdiags([-ones(N_img_pixels,1),ones(N_img_pixels,1)],[0,1],N_img_pixels,N_img_pixels);
last_row_indices_mask=false(N_img_pixels,1);last_row_indices_mask(rows:rows:N_img_pixels)=true;Dy(last_row_indices_mask,:)=0;
end