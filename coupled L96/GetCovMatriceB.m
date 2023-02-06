function [Bainv,Boinv,Ba,Bo,SD] = GetCovMatriceB(number_of_samples,h,nsteps,na,no,Fx,Fy,alph,gamma,l_SpCov_SOAR,...
        L_atmos, L_ocean,variance_atmos, variance_ocean)
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

if l_SpCov_SOAR == 1
    % SOAR types B:
    [Ba] = generateBforL96c(na,L_atmos,variance_atmos);
    [Bo] = generateBforL96c(no,L_ocean,variance_ocean);
    B = blkdiag(Ba,Bo);
    SD = [sqrt(variance_atmos)*ones(na,1); sqrt(variance_ocean)*ones(no,1)];
end
% Sample Covaraince from simulations
if l_SpCov_SOAR == 0
    % number_of_samples = floor(0.8*2*N);
    [B,C,s,SD] = sdcal_covgen_l96c([x0_init';y0_init'],h,nsteps,number_of_samples,na,no,Fx,Fy,alph,gamma);
    B(1:na,1:na) = variance_atmos*B(1:na,1:na);
    B(na+1:end,na+1:end) = variance_ocean*B(na+1:end,na+1:end);
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

