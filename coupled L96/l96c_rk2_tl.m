function [dz,dfdz_store,dfdYi_store] = l96c_rk2_tl(z_lin,dz0,h,nsteps,na,no,Fx,Fy,alph,gamma)
% Tangent Linear state forward marching using 2 stage Heun's RK2 method

% Initialise output variables
z = zeros(na+no,nsteps+1);
dz = zeros(na+no,nsteps+1);
Yi = zeros(na+no,nsteps);
dfdz_store = zeros(na+no,na+no,nsteps+1);
dfdYi_store = zeros(na+no,na+no,nsteps+1);
z = z_lin;
dz(:,1) = dz0;
%% rk2 tlm version: dz[n+1] = dz[n] + (h/2)*(k1 + k2); l1 = tlm(z[n])*dz[n], l2 = tlm[Yi]*dYi;
%% dYi = dz(:,i) + h*l1; Yi = z[n] + h*k1;

for i=1:nsteps    
    % rk2 update of the tlm:
    Yi = z(:,i) + h*L96_coupled_fn(z(:,i),na,Fx,Fy,alph,gamma);
    dfdz = tlm_l96c(z(:,i),no,na,alph,gamma);
    dfdz_store(:,:,i) = dfdz;
    dk1 = dfdz*dz(:,i);
    dYi = dz(:,i) + h*dk1;
    dfdYi = tlm_l96c(Yi,no,na,alph,gamma);
    dk2 = dfdYi*dYi;
    dz(:,i+1) = dz(:,i) + 0.5d0*h*(dk1+dk2);
    dfdYi_store(:,:,i) = dfdYi; % Jacob(Yi) is stored and output for adjoint computation later
end
end

