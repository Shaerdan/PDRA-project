function [g_hat_model] = l96c_rk2_adj(dfdz,dfdYi,g_hat,h,nsteps,na,no)
% Adjoint code for rk2 2 stage Heun's method, L96 coupled model.
%
% Input variables:
% Jacob1 = df(z)/dz; Jacob2 = df(Yi)/dYi; 
% g_hat is the observational gradient at one next future step (backward iteration)
%
% Output variables:
% g_hat_model is the descrete adjoint of rk2: 2 stage Heun's method, L96 coupled model. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g_hat_model = zeros(na+no,nsteps+1);
g_hat_model(:,nsteps+1) = g_hat;
% u = zeros(na+no,nsteps+1);
for i=nsteps:-1:1
    u2(:,i) = 0.5d0*h* dfdYi(:,:,i)'*(g_hat_model(:,i+1));    
    u1(:,i) = h* dfdz(:,:,i)'*(0.5d0*g_hat_model(:,i+1) + u2(:,i));
    g_hat_model(:,i) = g_hat_model(:,i+1) + u1(:,i)+ u2(:,i);
end  
end    
    
