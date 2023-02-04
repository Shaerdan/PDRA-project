function [f,g] = calcfg_ocean_l96c(dX0_a,z_b,innov,z_lin,H,invB,invR,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% compute the cost function for the atmospheric part
% x=z_lin(1:na,:);
% y=z_lin(na+1:end,:);

%% Calculate cost function
% Background term
dXb_o = z_b(na+1:end,1) - z_lin(na+1:end,1);
fa = 0.5 * (dX0_a-dXb_o)' * invB * (dX0_a-dXb_o);

z = z_lin;
dz0 = [zeros(na,1);dX0_a];
% [dz,M] = transmat_lor96c(z0,nsteps*h,h,na,Fx,Fy,alph,gamma);
[dz,dfdz,dfdYi] = l96c_rk2_tl(z,dz0,h,nsteps,na,no,Fx,Fy,alph,gamma);
% [dx,dy,dz,dw,dv] = molteni_rk2_tl(h,nsteps,beta,rho,sigma,alpha_b, ...
%     Omega,k,w_star,x,y,z,w,v,coupling_freq,dX0_o(1),dX0_o(2),dX0_o(3),0,0,alpha_b);
for i=1:nsteps
    fa = fa + 0.5 * ( ob_ix(na+1:end,i)'.*(innov(na+1:end,i)-H(na+1:end,na+1:end)*dz(na+1:end,i+1))'...
        *invR*(innov(na+1:end,i)-H(na+1:end,na+1:end)*dz(na+1:end,i+1)));
end
f = fa;

% initialising the observation gradient, compute observation gradient at the last step:
g_hat = zeros(na+no,nsteps+1);
g_hat(na+1:na+no,nsteps+1) = g_hat(na+1:na+no,nsteps+1) - invR*(ob_ix(na+1:na+no,nsteps).*...
    (innov(na+1:na+no,nsteps)-H(na+1:na+no,na+1:na+no)*dz(na+1:na+no,nsteps+1)));
% adjoint interation:
for i=nsteps:-1:2
%   One step of adjoint model  
    dfdz_adj = dfdz(:,:,i);
    dfdYi_adj = dfdYi(:,:,i);
    [g_hat_model] = l96c_rk2_adj(dfdz_adj,dfdYi_adj,g_hat(:,i+1),h,1,na,no);    
    g_hat(na+1:end,i) = g_hat_model(na+1:end,1) - invR*(ob_ix(na+1:end,i-1).*...
        (innov(na+1:end,i-1)-H(na+1:end,na+1:end)*dz(na+1:end,i)));
end
[g_hat_model] = l96c_rk2_adj(dfdz(:,:,1),dfdYi(:,:,1),g_hat(:,2),h,1,na,no);
g0 = g_hat_model(:,1);
% l96c_rk2_adjoint(Jacob_out1, Jacob_out2, Yi, z_final, dz_final, h, nsteps)
gJb = invB * (dX0_a-dXb_o);
g = gJb + g0(na+1:end);

end

