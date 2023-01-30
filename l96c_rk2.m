function [z] = l96c_rk2(z0,h,nsteps,na,no,Fx,Fy,alph,gamma)
% Integrate coupled model of Molteni et al. 1993 using 2nd order RK2
%
% Input variables:
% z0 is the initial state
%
% Output variables:
% z is the state stored evolved from t_1 to t_nsteps, all stored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise output variables
z = zeros(na+no,nsteps+1);
z(:,1) = z0;

%% rk2 version: z[n+1] = z[n] + (h/2)*(k1 + k2); k1 = f(z[n]), k2 = f[z[n]+h*k1];
%       
for i=1:nsteps
    k1 = f_model_l96c(z(:,i),no,na,alph,gamma,Fx,Fy);
    Yi = z(:,i) + h*k1;
    k2 = f_model_l96c(Yi,no,na,alph,gamma,Fx,Fy);
    z(:,i+1) = z(:,i) + 0.5d0*h*(k1+k2);
end
end
