function [f,g] = calcfg_atmos_l96c(dX0_a,z_b,innov,z_lin,H,invB,invR,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% compute the cost function for the atmospheric part
% x=z_lin(1:na,:);
% y=z_lin(na+1:end,:);

%% Calculate cost function
% Background term
dXb_a = z_b(1:na,1) - z_lin(1:na,1);
fa = 0.5 * (dX0_a-dXb_a)' * invB * (dX0_a-dXb_a);

z = z_lin;
dz0 = [dX0_a;zeros(no,1)];
% [dz,M] = transmat_lor96c(z0,nsteps*h,h,na,Fx,Fy,alph,gamma);
[dz,dfdz,dfdYi] = l96c_rk2_tl(z,dz0,h,nsteps,na,no,Fx,Fy,alph,gamma);
% [dx,dy,dz,dw,dv] = molteni_rk2_tl(h,nsteps,beta,rho,sigma,alpha_b, ...
%     Omega,k,w_star,x,y,z,w,v,coupling_freq,dX0_a(1),dX0_a(2),dX0_a(3),0,0,alpha_b);
for i=1:nsteps
    fa = fa + 0.5 * (ob_ix(1:na,i)'.*(innov(1:na,i)-H(1:na,1:na)*dz(1:na,i+1))'...
        *invR*((innov(1:na,i)-H(1:na,1:na)*dz(1:na,i+1))));
end
f = fa;

%% observation gradient at the last step:
g_hat = zeros(na+no,nsteps+1);
g_hat(1:na,nsteps+1) = g_hat(1:na,nsteps+1) - invR*(ob_ix(1:na,nsteps).*...
    (innov(1:na,nsteps)-H(1:na,1:na)*dz(1:na,nsteps+1)));
% adjoint interation:
for i=nsteps:-1:2
    dfdz_adj = dfdz(:,:,i);
    dfdYi_adj = dfdYi(:,:,i);
    %   One step of adjoint model
    [g_hat_model] = l96c_rk2_adj(dfdz_adj,dfdYi_adj,g_hat(:,i+1),h,1,na,no);
    g_hat(1:na,i) = g_hat_model(1:na,1) - invR*(ob_ix(1:na,i-1).*(innov(1:na,i-1)-H(1:na,1:na)*dz(1:na,i)));
end
[g_hat_model] = l96c_rk2_adj(dfdz(:,:,1),dfdYi(:,:,1),g_hat(:,2),h,1,na,no);
g0 = g_hat_model(:,1);
gJb = invB * (dX0_a-dXb_a);
g = gJb+ g0(1:na);

end

