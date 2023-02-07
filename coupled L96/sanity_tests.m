
clc; clear all; close all;
% tl test, adjoint test, gradient test:

plot_num =1 ;
h = 0.0125d0;
nsteps = 50;
na = 40;
no = 40;
n = na + no;
%% TLM Test and Adjoint Test:
%% test the tlm for atmosphere and ocean separately for the weakly coupled system, we do this
%% by setting gamma = 0:
Fx=15;
Fy=8;
alph=0.5;
gamma= 0.6;
number_of_samples = n; % full sample size for (likely) nonsingular B
l_SpCov_SOAR = 0;
L_atmos = 2; L_ocean = 4; variance_atmos = 0.1; variance_ocean = 0.1;
[Bainv,Boinv,Ba,Bo,B,SD] = GetCovMatriceB(number_of_samples,h,nsteps,na,no,Fx,Fy,alph,gamma,l_SpCov_SOAR,...
        L_atmos, L_ocean,variance_atmos, variance_ocean);
% Jacobian test:
% f1 = f_model_l96c(z0,no,na,alph,gamma,Fx,Fy);
z0 = rand(n,1);
dz0 = rand(n,1);
f1 = L96_coupled_fn(z0,na,Fx,Fy,alph,0);
dfdz = tlm_l96c(z0,na,no,alph,0);
for i=1:14
    a_jacobian_test(i) = 10^(-i);
    z0_perturbed = z0 + a_jacobian_test(i)*dz0;
    %     f2 = f_model_l96c(z0_perturbed,no,na,alph,gamma,Fx,Fy);
    f2 = L96_coupled_fn(z0_perturbed,na,Fx,Fy,alph,0);
    jac_test(i) = norm(f2 - f1)/norm(a_jacobian_test(i)*dfdz*dz0);
end
figure(plot_num)
loglog(a_jacobian_test,abs(jac_test-1),'DisplayName','Jacobian test for weakly coupled');
xlabel('perturbation scaling')
ylabel('$\frac{||f(z_0+\gamma dz_0)-f(z_)||}{||\gamma*\frac{df}{dz}|_{z_{0}}*dz_0||} - 1$','interpreter','latex')
legend show
plot_num = plot_num + 1;

z0 = rand(n,1);
dz0 = rand(n,1);
% calculate zk = M_nonlinear(z0), dz_k = M_linear(z0,dz0)
z_lin = l96c_rk2(z0,h,nsteps,na,no,Fx,Fy,alph,0);
z=z_lin;
% note that z stores the state at all time steps, for the test we only need
% the last step of z:
zk = z(:,end);
% dz0 = z(:,2) - z(:,1);
[dz,Jacob_out1,Jacob_out2] = l96c_rk2_tl(z_lin,dz0,h,nsteps,na,no,Fx,Fy,alph,0);
% similar to the case of z, we only need dz(:,end):
dzk = dz(:,end);
dzk_norm = norm(dzk);
for i=1:14
    a_tl(i) = 10^(-i);
    z0_perturbed = z0 + a_tl(i)*dz0;
    z_perturbed = l96c_rk2(z0_perturbed,h,nsteps,na,no,Fx,Fy,alph,0);
    zk_perturbed = z_perturbed(:,end);
    tl_test_weakly(i) = norm(zk_perturbed - zk - a_tl(i)*dzk)/(a_tl(i)*dzk_norm);
end
figure(plot_num)
loglog(a_tl,tl_test_weakly,'k-*','DisplayName','tlm test')
xlabel('perturbation scaling')
ylabel('error')
legend show
plot_num = plot_num +1;

% Adjoint test:
[g_hat_model] = l96c_rk2_adj(Jacob_out1,Jacob_out2,dzk,h,nsteps,na,no);
lhs = dzk'*dzk;
rhs = dz0'*g_hat_model(:,1);
adj_test_weakly = lhs - rhs;
format long
disp(strcat('adjoint test, ',' error =  ',num2str(adj_test_weakly)))

%% Gradient Test:
% for gradient test, set gamma back to the nonzero value for weakly coupled
% system.
gamma = 0.6;

%% Gradient Tests For calcfg_atmos_l96c and calcfg_ocean_l96c:
H = eye(n,n);
% Bainv = rand(na,na); Bainv = Bainv'*Bainv;
% Boinv = rand(no,no); Boinv = Boinv'*Boinv;
Rainv = rand(na,na); Rainv = Rainv'*Rainv;
Roinv = rand(no,no); Roinv = Roinv'*Roinv;
% ob_ix = ones(n,nsteps);
n_cycles_per_smoother = 4; % number of short window cycles per smoother cycle
n_ob_pattern_repeats = 4; % Total length of the run in terms of ob_pattern_repeat_freq
ob_pattern_repeat_freq = 1; % The pattern of observation times repeats after
assim_steps =n_cycles_per_smoother*nsteps;
z_lin = l96c_rk2(z0,h,assim_steps,na,no,Fx,Fy,alph,gamma);
z=z_lin;

ob_ix=zeros(n,assim_steps);
obs_noise=randn(n,assim_steps);
u_ob=zeros(n,assim_steps);
i_part_of_ob_pattern = 1;
x_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
y_ob = (ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps:1:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
x_ob_local = x_ob((x_ob > (i_part_of_ob_pattern-1)*assim_steps) & (x_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
y_ob_local = y_ob((y_ob > (i_part_of_ob_pattern-1)*assim_steps) & (y_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
var_ob = [1 1];
% obs_noise = randn(na+no,1);
for i=x_ob_local
    ob_ix(i) = 1;
    u_ob(1:na,i) = z_lin(1:na,i+1) + sqrt(var_ob(1))*obs_noise(1:na,i);
end
for i=y_ob_local
    ob_ix(i) = 1;
    u_ob(na+1:n,i) =  z_lin(na+1:n,i+1) + sqrt(var_ob(2))*obs_noise(na+1:n,i);
end
%
%
dX0_a = rand(na,1); dX0_o = rand(no,1); 
ub = zeros(n,1);innov=rand(n,assim_steps); u_lin = z_lin;
u_lin(:,1) = ub; 
[fa1,ga1] = calcfg_atmos_l96c(dX0_a,ub,innov,...
    u_lin,H,Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix);
hg_a = ga1/norm(ga1,2);
[fo1,go1] = calcfg_ocean_l96c(dX0_o,ub,innov,...
    u_lin,H,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix);
hg_o = go1/norm(go1,2);
for itest = 1:14
    a(itest) = 10^(-itest);
    [fa2,ga2] = calcfg_atmos_l96c(dX0_a+a(itest)*hg_a,ub,innov,...
        u_lin,H,Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix);
    Phi_a_atmos(itest) = (fa2 - fa1)/(a(itest)*hg_a'*ga1);
    [fo2,go2] = calcfg_ocean_l96c(dX0_o+a(itest)*hg_o,ub,innov,...
        u_lin,H,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix);
    Phi_a_ocean(itest) = (fo2 - fo1)/(a(itest)*hg_o'*go1);
end
figure(150)
loglog(a,abs(Phi_a_atmos-1),'k','DisplayName','Gradient Test for Atmospheric Costs and Gradient')
xlabel('Perturbeation Scaling')
ylabel('$|\frac{||J(dXa+\gamma*h) - J(dXa)||}{||\gamma*h*dJ||}-1|$','interpreter','latex')
legend show

figure(151)
loglog(a,abs(Phi_a_ocean-1),'k','DisplayName','Gradient Test for Oceanic Costs and Gradient')
xlabel('Perturbeation Scaling')
ylabel('$|\frac{||J(dXa+\gamma*h) - J(dXa)||}{||\gamma*h*dJ||}-1|$','interpreter','latex')
legend show

% Test the nonlinear model:
nsteps_long = 50;
h_model_test = 0.0125d0;
xvals=1:na; % atmosphere grid indeces
yvals=1:no; % ocean grid indeces
x0_init=sin(xvals/(na-1)*2*pi);
y0_init=cos(5*yvals/(no-1)*2*pi);

z0_model_test = [x0_init';y0_init']+50*rand(n,1);
[z_model_test] = l96c_rk2(z0_model_test,h_model_test,nsteps_long,na,no,Fx,Fy,alph,gamma);

%     case 'fully coupled'
%         %% here we also can test the tlm for the fully coupled model(for future applicaitons):
%         Fx=15;
%         Fy=8;
%         alph=0.5;
%         gamma= 0.6;
%
%
%         % Jacobian test:
%         % f1 = f_model_l96c(z0,no,na,alph,gamma,Fx,Fy);
%         z0 = rand(n,1);
%         dz0 = rand(n,1);
%         f1 = L96_coupled_fn(z0,na,Fx,Fy,alph,gamma);
%         dfdz = tlm_l96c_fullcouple(z0,na,no,alph,gamma);
%         for i=1:14
%             a_jacobian_test(i) = 10^(-i);
%             z0_perturbed = z0 + a_jacobian_test(i)*dz0;
%             %     f2 = f_model_l96c(z0_perturbed,no,na,alph,gamma,Fx,Fy);
%             f2 = L96_coupled_fn(z0_perturbed,na,Fx,Fy,alph,gamma);
%             jac_test(i) = norm(f2 - f1)/norm(a_jacobian_test(i)*dfdz*dz0);
%         end
%         figure(plot_num)
%         loglog(a_jacobian_test,abs(jac_test-1),'DisplayName','fully coupled');
%         xlabel('perturbation scaling')
%         ylabel('$\frac{||f(z_0+\gamma dz_0)-f(z_)||}{||\gamma*\frac{df}{dz}|_{z_{0}}*dz_0||} - 1$','interpreter','latex')
%         legend show
%         plot_num = plot_num + 1;
%
%
%
%         z0 = rand(n,1);
%         dz0 = rand(n,1);
%         % calculate zk = M_nonlinear(z0), dz_k = M_linear(z0,dz0)
%         z_lin = l96c_rk2(z0,h,nsteps,na,no,Fx,Fy,alph,gamma);
%         z = z_lin;
%         % note that z stores the state at all time steps, for the test we only need
%         % the last step of z:
%         clear Jacob_out1 Jacob_out2
%         zk = z(:,end);
%         [dz,Jacob_out1,Jacob_out2] = l96c_rk2_tl_fullcouple(z_lin,dz0,h,nsteps,na,no,Fx,Fy,alph,gamma);
%         % similar to the case of z, we only need dz(:,end):
%         dzk = dz(:,end);
%         dzk_norm = norm(dzk);
%         for i=1:12
%             a_tl_full(i) = 10^(-i);
%             z0_perturbed = z0 + a_tl_full(i)*dz0;
%             z_purterbed = l96c_rk2(z0_perturbed,h,nsteps,na,no,Fx,Fy,alph,gamma);
%             zk_purterbed = z_purterbed(:,end);
%             tl_test_fully(i) = norm(zk_purterbed - zk - a_tl_full(i)*dzk)/(a_tl_full(i)*dzk_norm);
%         end
%         [g_hat_model] = l96c_rk2_adj(Jacob_out1,Jacob_out2,dzk,h,nsteps,na,no);
%         lhs = dzk'*dzk;
%         rhs = dz0'*g_hat_model(:,1);
%         adj_test_fully = lhs - rhs;
%         disp(strcat('adjoint test',num2str(adj_test_fully)))
%         figure(plot_num)
%         loglog(a_tl_full,tl_test_fully,'k-*','DisplayName','tlm test for fully coupled model')
%         xlabel('perturbation scaling')
%         ylabel('error')
%         legend show
%         plot_num = plot_num +1;
%         disp(strcat('Adjoint test results ',num2str(adj_test_weakly),',',num2str(adj_test_fully)))
% end
