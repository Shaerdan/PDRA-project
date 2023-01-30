function [B,C,s,SD] = sdcal_covgen_l96c(z0,h,nsteps,number_of_samples,na,no,Fx,Fy,alph,gamma)
%
% (1) Calculate natural variability of the variables, to be used as normalisation factors for the error norms
% (2) Generate covariance matrix using Canadian Quick Covariances method
%
% Total number of steps
total_steps = nsteps * number_of_samples;

[z] = l96c_rk2(z0,h,total_steps,na,no,Fx,Fy,alph,gamma);

% Standard deviation of the variables
SD = std(z');

% Generate samples for covariance calculation
j = 1;
for i=nsteps+1:nsteps:total_steps+1
    dz(:,j) = z(:,i) - z(:,i-nsteps);
    j = j+1;
end
X = dz';
B = cov(X);
[C,s]=corrcov(B);