function [Bainv,Boinv,Ba,Bo,SD] = GetCovMatriceB(number_of_samples,h,nsteps,na,no,Fx,Fy,alph,gamma)
%UNTITLED Summary of this function goes here
%   formulate covariance matrices
%   R matrix need some updating

xvals=1:na; % atmosphere grid indeces
yvals=1:no; % ocean grid indeces
x0_init=sin(xvals/(na-1)*2*pi);
y0_init=cos(5*yvals/(no-1)*2*pi);
% fun experiments with the model for N<=4:
% [z_chk] = l96c_rk2([x0_init';y0_init'],h,nsteps,na,no,Fx,Fy,alph,gamma);
% figure(900)
% plot3(z_chk(1+na,:),z_chk(2+na,:),z_chk(3+na,:),'r-'); hold on; ...
%     plot3(z_chk(1,:),z_chk(2,:),z_chk(3,:),'k-');

% y_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);

%% formulate background&observation error cov matrices, and H matrix:

B_method = 1;
if B_method == 0
    % SOAR types B:
    L_atmos = 2; L_ocean = 4; % make these input variable
    variance_atmos = 1; variance_ocean = 1; % make these input variable
    [Ba] = generateBforL96c(N,L_atmos,variance_atmos);
    [Bo] = generateBforL96c(N,L_ocean,variance_ocean);
    B = blkdiag(Ba,Bo);
end
% Sample Covaraince from simulations
if B_method == 1
    % number_of_samples = floor(0.8*2*N);
    [B,C,s,SD] = sdcal_covgen_l96c([x0_init';y0_init'],h,nsteps,number_of_samples,na,no,Fx,Fy,alph,gamma);
    Ba = B(1:na,1:na);
    Bo = B(na+1:end,na+1:end);
    figure(300)
    imagesc(B);
end
Bainv = inv(Ba);
Boinv = inv(Bo);
% Bainv = eye(na,na);
% Boinv = eye(no,no);

end

