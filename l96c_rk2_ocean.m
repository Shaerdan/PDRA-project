function [y] = l96c_rk2_ocean(x,y0,h,nsteps,no,Fy,alph,gamma)
% Integrate coupled model of Molteni et al. 1993 using 2nd order RK2
%
% Input variables:

%
% Output variables:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise output variables
y = zeros(no,nsteps+1);
y(:,1) = y0;

%% rk2 version: z[n+1] = z[n] + (h/2)*(k1 + k2); k1 = f(z[n]), k2 = f[z[n]+h*k1];
%       
for i=1:nsteps
    x_c(:,i) = x(:,i);
    k1 = f_model_l96c_ocean(y(:,i),x_c,no,alph,gamma,Fy);
    Yi = y(:,i) + h*k1;
    k2 = f_model_l96c_ocean(Yi,x_c,no,alph,gamma,Fy);
    y(:,i+1) = y(:,i) + 0.5d0*h*(k1+k2);
end
end
