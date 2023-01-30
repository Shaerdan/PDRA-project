function [dz,Jacob_out1,Jacob_out2] = l96c_rk2_tl(z0,dz0,h,nsteps,na,no,Fx,Fy,alph,gamma)
% Tangent Linear state forward marching using 2 stage Heun's RK2 method

% Initialise output variables
z = zeros(na+no,nsteps+1);
dz = zeros(na+no,nsteps+1);
Yi = zeros(na+no,nsteps);
Jacob_out1 = zeros(na+no,na+no,nsteps+1);
Jacob_out2 = Jacob_out1;
z(:,1) = z0;
dz(:,1) = dz0;

%% rk2 tlm version: dz[n+1] = dz[n] + (h/2)*(k1 + k2); l1 = tlm(z[n])*dz[n], l2 = tlm[Yi]*dYi;
%% dYi = dz(:,i) + h*l1; Yi = z[n] + h*k1;
for i=1:nsteps
% rk2 update of the state:
    k1 = f_model_l96c(z(:,i),no,na,alph,gamma,Fx,Fy);
    Yi(:,i) = z(:,i) + h*k1; % Yi is stored and output for adjoint computation later
    k2 = f_model_l96c(Yi,no,na,alph,gamma,Fx,Fy);
    z(:,i+1) = z(:,i) + 0.5d0*h*(k1+k2);
% rk2 update of the tlm:
    Jacob1 = tlm_l96c(z(:,i),no,na,alph,gamma);
    Jacob_out1(:,:,i) = Jacob1;    
    l1 = Jacob1*dz(:,i);
    dYi = dz(:,i) + h*l1;
    Jacob2 = tlm_l96c(Yi,no,na,alph,gamma); 
    l2 = Jacob2*dYi;
    dz(:,i+1) = dz(:,i) + 0.5d0*h*(l1+l2);
    Jacob_out2(:,:,i) = Jacob2; % Jacob(Yi) is stored and output for adjoint computation later
end
end

