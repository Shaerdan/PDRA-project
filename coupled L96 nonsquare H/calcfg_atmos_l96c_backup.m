function [f,g] = calcfg_atmos_l96c(dX0_a,u_b,innov,u_lin,H,invB,invR,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% compute the cost function for the atmospheric part
% x=u_lin(1:na,:);
% y=u_lin(na+1:end,:);

%% Calculate cost function
% Background term
dXb_a = u_b(1:na,1) - u_lin(1:na,1);
fa = 0.5 * (dX0_a-dXb_a)' * invB * (dX0_a-dXb_a);

z0 = u_lin(:,1);
dz0 = [dX0_a;zeros(no,1)];
% [dz,M] = transmat_lor96c(z0,nsteps*h,h,na,Fx,Fy,alph,gamma);
[dz,Jacob1,Jacob2] = l96c_rk2_tl(z0,dz0,h,nsteps,na,no,Fx,Fy,alph,0);
% [dx,dy,dz,dw,dv] = molteni_rk2_tl(h,nsteps,beta,rho,sigma,alpha_b, ...
%     Omega,k,w_star,x,y,z,w,v,coupling_freq,dX0_a(1),dX0_a(2),dX0_a(3),0,0,alpha_b);
for i=1:nsteps
    fa = fa + 0.5 * ( ob_ix(1:na,i)'.*(innov(1:na,i)-H(1:na,1:na)*dz(1:na,i+1))'...
        *invR*(innov(1:na,i)-H(1:na,1:na)*dz(1:na,i+1)));
end
f = fa;
g_hat = zeros(na+no,nsteps+1);
% observation gradient at the last step:
g_hat(1:na,nsteps+1) = g_hat(1:na,nsteps+1) - invR*(ob_ix(1:na,nsteps).*...
    (innov(1:na,nsteps)-H(1:na,1:na)*dz(1:na,nsteps+1)));
Jacob1_adj = zeros(no+na,no+na); Jacob2_adj = Jacob1_adj;

[g_hat_model] = l96c_rk2_adj(Jacob1_adj,Jacob2_adj,dz(:,nsteps+1),h,1,na,no);    

for i=nsteps:-1:2
%     Jacob1_adj(1:na,1:na) = Jacob1(1:na,1:na,i);
%     Jacob2_adj(1:na,1:na) = Jacob2(1:na,1:na,i);
%   One step of adjoint model  
%     [g_hat_model] = l96c_rk2_adj(Jacob1_adj,Jacob2_adj,g_hat(:,i+1),h,1,na,no);    
%     [g_hat_model] = l96c_rk2_adjoint(Jacob1, Jacob2, z_final, dz_final, h, 1); 
%     g_hat(1:na,i) = g_hat_model(1:na,1) - invR*(ob_ix(1:na,i-1).*(innov(1:na,i-1)-H(1:na,1:na)*dz(1:na,i)));
g_obs(:,i) = invR*(ob_ix(1:na,i-1).*(innov(1:na,i-1)-H(1:na,1:na)*dz(1:na,i)))*g_hat_model(:,i);
end
g0 = sum(g_obs);
%     [g_hat_model] = l96c_rk2_adj(Jacob1_adj,Jacob2_adj,g_hat(:,2),h,1,na,no);        
%     g0 = g_hat_model(:,1);
% l96c_rk2_adjoint(Jacob_out1, Jacob_out2, Yi, z_final, dz_final, h, nsteps)
gJb = invB * (dX0_a-dXb_a);
g = gJb+ g0(1:na);

end

