clc; clear all; close all;
tolerance = 1.0d-10;   % For solver
max_iterations = 200; % For solver
% Grad_Test = 1: turn on gradient tests for the adjoint code; = 0 turn off.
Grad_Test = 0;
%% rk2 solver parameters and model parameters for the l96 coupled model
nsteps = 50;
h=0.00125;
Fx=15;
Fy=8;
alph=0.5;
gamma= 0.6;
N = 40;
na = N; no = N; ntotal = na + no;
% loop controls:
outer_loops=2;     % number of outerloop for weakly coupled standard 4dvar
n_cycles_per_smoother = 4; % number of short window cycles per smoother cycle
n_ob_pattern_repeats = 10; % Total length of the run in terms of ob_pattern_repeat_freq
ob_pattern_repeat_freq = 1; % The pattern of observation times repeats after
% method control:
assim_scheme = 5;  % 5 for smoother method
l_fgat_s5 = 0;     % 0 = 4DVar for smoother step; 1 = 3DFGAT for smoother step
% data control:
l_newbg_xb = 1;
l_newobs = 1;
% plot control:
l_plot_convergence = 0;
l_plot_state = 0;
l_plot_error_norm = 1; 
l_plot_avg_error_norm = 0;
l_plot_trajectories = 1;
%% Smoother setup
s5_B_scaling = 1;
s5_smoother_loops = 2;  % Number of outer loops for smoother step only
s5_iterations = 2;
l_lin_s5 = 1;       % 0 = Take the analysis trajectory as both background and the first linearisation state;
l_integration_coupled_s5 = 1;   % At the last outer loop of the last iteration of the smoother step (which produces
% the final analysis), whether to integrate the smoother analysis at initial time using
% the coupled model; if yes, the atmospheric trajectory would be reset to that of the
% original atmospheric analyses at the beginning of each standard assimilation window
% before the coupled model integration is continued (this capability hasn't been designed
% for cases where IAU is used for the smoother step)
update_method = 1;
reset_atmosphere = 1;
Increment_Scaling = 1.0d0;
%% Setting up parameters for the assimilation

% this number of (smoother) assimilation cycles
assim_steps = nsteps*n_cycles_per_smoother;

data_bgx_out='data_bgBx.mat';  % Output file for xb
data_obs_out='data_obs.mat';   % Output file or x_obs


xvals=1:na; % atmosphere grid indeces
yvals=1:no; % ocean grid indeces
x0_init=sin(xvals/(na-1)*2*pi);
y0_init=cos(5*yvals/(no-1)*2*pi);
% fun experiments with the model for N<=4:
[z_chk] = l96c_rk2([x0_init';y0_init'],h,10*assim_steps,na,no,Fx,Fy,alph,gamma);
figure(900)
plot3(z_chk(1+na,:),z_chk(2+na,:),z_chk(3+na,:),'r-'); hold on; ...
plot3(z_chk(1,:),z_chk(2,:),z_chk(3,:),'k-'); 

x_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
y_ob = (ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps:1:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
% y_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
figure(400)
imagesc(x_ob)

%% formulate background&observation error cov matrices, and H matrix:
L_atmos = 2; L_ocean = 4;
variance_atmos = 0.1; variance_ocean = 0.1;
B_method = 1;
if B_method == 0
% SOAR types B:
[Bc_atmos] = generateBforL96c(N,L_atmos,variance_atmos);
[Bc_ocean] = generateBforL96c(N,L_ocean,variance_ocean);
B = blkdiag(Bc_atmos,Bc_ocean);
end
% Sample Covaraince from simulations
if B_method == 1
% number_of_samples = floor(0.8*2*N);
number_of_samples = 80;
[B,C,s,SD] = sdcal_covgen_l96c([x0_init';y0_init'],h,nsteps,number_of_samples,na,no,Fx,Fy,alph,gamma);
Bc_atmos = B(1:na,1:na);
Bc_ocean = B(na+1:end,na+1:end);
figure(300)
imagesc(B);
end
% plot(Bc_atmos(20,:))
% plot(Bc_atmos(20,:))
Bainv = inv(Bc_atmos);
Boinv = inv(Bc_ocean);
% Binv = inv(B);
var_ob = [1.0^-4, 1.0^-4];
R_atmos = var_ob(1)*eye(na,na); R_ocean = var_ob(2)*eye(no,no);
R = blkdiag(R_atmos,R_ocean);
Rinv = inv(R);
Rainv = inv(R_atmos);
Roinv = inv(R_ocean);
H = blkdiag(var_ob(1)*eye(na,na),var_ob(1)*eye(no,no));
for i_ob_pattern_repeats = 1:n_ob_pattern_repeats
    for i_part_of_ob_pattern = 1:ob_pattern_repeat_freq
        %% Start identical-twin set-up
        
        ub_plot=zeros(ntotal,n_cycles_per_smoother,nsteps+1);
        
        % Generate truth
        if (i_ob_pattern_repeats == 1 && i_part_of_ob_pattern == 1)
            x0_t = x0_init;
            y0_t = y0_init;
        else
            x0_t = x(:,assim_steps+1)';
            y0_t = y(:,assim_steps+1)';
        end
        z0_t = [x0_t,y0_t];
        %         [z,M] = transmat_lor96c(z0_t,assim_steps*h,h,na,Fx,Fy,alph,gamma);
        [z] = l96c_rk2(z0_t,h,assim_steps,na,no,Fx,Fy,alph,gamma);
        %         [dztest] = l96c_rk2_tl(z0_t,dz0,h,assim_steps,na,no,Fx,Fy,alph,gamma);
        
        u = z;
        x = z(1:na,:);
        y = z(na+1:end,:);
        u_t=[x0_t,y0_t]';
        X_t= u_t;
        
        % Generate background
        if (i_ob_pattern_repeats == 1 && i_part_of_ob_pattern == 1)
            if l_newbg_xb
                noise = randn(ntotal,1);
                u_b = u_t + sqrtm(B) * noise(1:ntotal);
                save(data_bgx_out,'u_b')
            else
                load(data_bgx_in,'u_b')
            end
        else
            u_b = ua_plot(:,end,end);
        end
        % initialising ua_plot
        ua_plot=zeros(ntotal,n_cycles_per_smoother,nsteps+1);
        
        if l_newobs
            % Generate obs
            ob_ix=zeros(ntotal,assim_steps);
            obs_noise=randn(ntotal,assim_steps);
            u_ob=zeros(ntotal,assim_steps);
            
            x_ob_local = x_ob((x_ob > (i_part_of_ob_pattern-1)*assim_steps) & (x_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
            y_ob_local = y_ob((y_ob > (i_part_of_ob_pattern-1)*assim_steps) & (y_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
            
            for i=x_ob_local
                ob_ix(1:na,i) = 1;
                u_ob(1:na,i) = z(1:na,i+1) + sqrt(var_ob(1))*obs_noise(1:na,i);
            end
            for i=y_ob_local
                ob_ix(na+1:ntotal,i) = 1;
                u_ob(na+1:ntotal,i) =  z(na+1:ntotal,i+1) + sqrt(var_ob(2))*obs_noise(na+1:ntotal,i);
            end
            
            save(data_obs_out,'u_ob','ob_ix')
            
            
        else
            disp('add input from saved data');
            %             load(data_obs_in)
        end
        %
        
        
        %% Outer loops for smoother
        
        for i_smooth_iteration=1:s5_iterations
            for i_cycles = 1:n_cycles_per_smoother
                if i_smooth_iteration == 1
                    %% Background forecast
                    disp('********* START OF NEW CYCLE *********')
                    X_b=u_b;
                    X_t=u_t;
                    Errormat_Tot_Atm_Oce_anal= zeros(s5_iterations, outer_loops, 3);
                    % h,nsteps+fcsteps,beta,rho,sigma,alpha_b, ...
                    %    Omega,k,w_star,u_b(1),u_b(2),u_b(3),u_b(4),u_b(5),coupling_freq
                    %                     [zb_f,M] = transmat_lor96c(u_b,nsteps*h,h,na,Fx,Fy,alph,gamma);
                    [zb_f] = l96c_rk2(u_b,h,nsteps,na,no,Fx,Fy,alph,gamma);
                    ub_f=zb_f;
                    xb_f = ub_f(1:na,:);
                    yb_f = ub_f(na+1:ntotal,:);
                    ub_plot(:,i_cycles,:) = ub_f;
                end
                
                %% Perform assimilation of the weakly coupled standard 4DVAR:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i_count = 1:outer_loops
                    % Calculate innovations using coupled model
                    if (i_smooth_iteration == 1 && i_count == 1)
                        % Use background trajectory
                        u0 = u_b;
                        x_lin = xb_f(:,1:nsteps+1);
                        y_lin = yb_f(:,1:nsteps+1);
                        
                    else
                        % Calculate nonlinear trajectory with coupled model
                        %                         [zb_f,M] = transmat_lor96c(u0,nsteps*h,h,na,Fx,Fy,alph,gamma);
                        [zb_f] = l96c_rk2(u0,h,nsteps,na,no,Fx,Fy,alph,gamma);
                        
                        %                         [x_lin,y_lin,z_lin,w_lin,v_lin] = molteni_rk2(h,nsteps,beta,rho,sigma,alpha_b, ...
                        %                             Omega,k,w_star,u0(1),u0(2),u0(3),u0(4),u0(5),coupling_freq);
                    end
                    u_lin = zb_f;
                    %                     u_lin = [x_lin,y_lin,z_lin,w_lin,v_lin]';
                    innov = u_ob(:,(i_cycles-1)*nsteps+1:i_cycles*nsteps) - H*u_lin(:,2:nsteps+1); % Size 5 x nsteps. Assume no ob at time zero.
                    % Includes where there are no obs, but can use ob_ix
                    % Perform separate minimisations
                    
                    %% weakly coupled standard 4dvar, minimisation routine:
                    % Atmosphere
                    dX0_a=zeros(na,1);
                    [dXa_anal, JXaInner, dJXaInner, ita, rel_grad_a] = minimize_mod_crit_NKN(dX0_a,'calcfg_atmos_l96c',max_iterations,tolerance, ...
                        ub_plot(:,i_cycles,1),innov,u_lin,H,...
                        Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                    % Ocean
                    dX0_o=zeros(no,1);
                    [dXo_anal, JXoInner, dJXoInner, ito, rel_grad_o] = minimize_mod_crit_NKN(dX0_o,'calcfg_ocean_l96c',max_iterations,tolerance, ...
                        ub_plot(:,i_cycles,1),innov,u_lin,H,...
                        Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                    
                    %%%%%%%%%%%%
                    format short
                    IT_Outer_ITa_Ito = [i_count, ita, ito];
                    format longe
                    rel_grad_a;
                    rel_grad_o;
                    Cost_RG_Atmos = [JXaInner, rel_grad_a];
                    Cost_RG_Ocean = [JXoInner, rel_grad_o];
                    %%%%%%%%%%%
                    
                    % Update fields
                    delta_u0 = [dXa_anal' dXo_anal']';
                    u0 = u0 + delta_u0;
                    
                    %%%%%%%%%%%
                    format long
                    X_anal = u0;
                    Errormat_Tot_Atm_Oce_anal(i_smooth_iteration, i_count, :)=...
                        [norm(X_anal(1:na)-X_t(1:na)),norm(X_anal(1:na)-X_t(1:na)),norm(X_anal(na+1:ntotal)-X_t(na+1:ntotal))];
                    %%%%%%%%%%
                    
                    if i_count == 1
                        JXa  = JXaInner;
                        dJXa = dJXaInner;
                        JXo  = JXoInner;
                        dJXo = dJXoInner;
                    else
                        JXa  = [JXa' JXaInner']';
                        dJXa = [dJXa dJXaInner];
                        JXo  = [JXo' JXoInner']';
                        dJXo = [dJXo dJXoInner];
                    end
                    if Grad_Test == 1
                        %% Gradient Tests For calcfg_atmos_l96c and calcfg_ocean_l96c:
                        dX0_a = randn(na,1); dX0_o = randn(no,1);
                        [fa1,ga1] = calcfg_atmos_l96c(dX0_a,ub_plot(:,i_cycles,1),innov,...
                            u_lin,H,Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                        hg_a = ga1/norm(ga1,2);
                        [fo1,go1] = calcfg_ocean_l96c(dX0_o,ub_plot(:,i_cycles,1),innov,...
                            u_lin,H,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                        hg_o = go1/norm(go1,2);
                        for itest = 1:14
                            a(itest) = 10^(-itest);
                            [fa2,ga2] = calcfg_atmos_l96c(dX0_a+a(itest)*hg_a,ub_plot(:,i_cycles,1),innov,...
                                u_lin,H,Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                            Phi_a_atmos(itest) = (fa2 - fa1)/(a(itest)*hg_a'*ga1);
                            [fo2,go2] = calcfg_ocean_l96c(dX0_o+a(itest)*hg_o,ub_plot(:,i_cycles,1),innov,...
                                u_lin,H,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                            Phi_a_ocean(itest) = (fo2 - fo1)/(a(itest)*hg_o'*go1);
                        end
                        figure(150)
                        loglog(a,abs(Phi_a_atmos-1),'b',a,abs(Phi_a_ocean-1),'k')
                        
                    end
                end % end outerloop
                if (i_cycles <= 4)
                    dXa_anal_cycle(:,i_cycles) = u0(1:na) - X_b(1:na);
                    %                                 Count_dXa = Count_dXa + 1;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %% Display error norms and store background errors for diagnostics
                disp('***Final Results***')
                format long
                Xb_Xt_Xanal =[X_b(1), X_t(1), X_anal(1)]
                %                 Bg_err(((i_ob_pattern_repeats-1) * ob_pattern_repeat_freq + (i_part_of_ob_pattern-1)) * n_cycles_per_smoother + i_cycles, :) = (X_b(1:5) - X_t(1:5))';
                %                 disp('Errors: Total, Atmos, Ocean')
                %                 error_b=[norm(rdivide(X_b(1:5)-X_t(1:5),SD')),norm(rdivide(X_b(1:3)-X_t(1:3),SD(1:3)')),norm(rdivide(X_b(4:5)-X_t(4:5),SD(4:5)'))]
                %                 error_anal=[norm(rdivide(X_anal(1:5)-X_t(1:5),SD')),norm(rdivide(X_anal(1:3)-X_t(1:3),SD(1:3)')),norm(rdivide(X_anal(4:5)-X_t(4:5),SD(4:5)'))]
                %%%%%%%%%%%%
                
                if l_plot_convergence
                    %                     %% Plot convergence
                    str_plot = ['Cycle ' num2str(i_cycles)];
                    %
                    plot_num = 100 + i_part_of_ob_pattern;
                    plot1=figure(plot_num);
                    %                      if (assim_scheme == 3 || assim_scheme == 4 || assim_scheme == 5 || assim_scheme == 6)
                    ja = max(size(JXa));
                    xvals=(0:ja-1);
                    subplot(n_cycles_per_smoother,4,(i_cycles-1)*4+1)
                    semilogy(xvals,JXa)
                    %    title('Atmosphere: Convergence of cost function')
                    title('Atmosphere: Cost')
                    xlabel('Iteration')
                    ylabel({str_plot; 'Cost function'} )
                    subplot(n_cycles_per_smoother,4,(i_cycles-1)*4+2)
                    grad=zeros(ja,1);
                    for i = 1:ja
                        g = dJXa(:,i);
                        grad(i) = norm(g);
                    end
                    semilogy(xvals,grad)
                    %    title('Atmosphere: Convergence of gradient')
                    title('Atmosphere: Grad')
                    xlabel('Iteration')
                    ylabel('Norm of gradient')
                    %
                    jo = max(size(JXo));
                    xvals=(0:jo-1);
                    subplot(n_cycles_per_smoother,4,(i_cycles-1)*4+3)
                    semilogy(xvals,JXo)
                    %    title('Ocean: Convergence of cost function')
                    title('Ocean: Cost')
                    xlabel('Iteration')
                    ylabel('Cost function')
                    subplot(n_cycles_per_smoother,4,(i_cycles-1)*4+4)
                    grad=zeros(jo,1);
                    for i = 1:jo
                        g = dJXo(:,i);
                        grad(i) = norm(g);
                    end
                    semilogy(xvals,grad)
                    %    title('Ocean: Convergence of gradient')
                    title('Ocean: Grad')
                    xlabel('Iteration')
                    ylabel('Norm of gradient')
                end
                
                %% Run forecast
                xa=X_anal(1:na);
                ya=X_anal(na+1:ntotal);
                
                za_f_ini = [xa;ya];
                [za_f] = l96c_rk2(za_f_ini,h,nsteps,na,no,Fx,Fy,alph,gamma);
                
                
                %                     [xa_f,ya_f,za_f,wa_f,va_f] = molteni_rk2(h,nsteps+fcsteps,beta,rho,sigma,alpha_anal, ...
                %                         Omega,k,w_star,xa,ya,za,wa,va,coupling_freq);
                
                ua_f=za_f;
                ua_plot(:,i_cycles,:) = ua_f;
                
                if i_smooth_iteration == 1
                    %% Set background and truth for next cycle
                    u_b = ua_f(:,nsteps+1);
                    u_t = z(:,i_cycles*nsteps+1);
                end
                if l_plot_state == 1
                    str_plot = ['Cycle ' num2str(i_cycles)];
                    plot_num2 = 10 + i_part_of_ob_pattern;
                    plot2=figure(plot_num2);
                    subplot(n_cycles_per_smoother/2,n_cycles_per_smoother/2,i_cycles)
                    plot(u_b,'k-','DisplayName','Analysis Forecast'); hold on;...
                        plot(u_t,'b-','DisplayName','True State'); hold on;
                    xlabel('Iteration')
                    ylabel({str_plot; 'soln'} )
                    legend show
                end
                
            end % i_cycles
            
            %% Smoother
            if (assim_scheme == 5)
                ua_plot_2 = reshape(permute(ua_plot(:,:,1:nsteps),[1,3,2]),[na+no,assim_steps]);
                ua_plot_2 = [ua_plot_2 ua_plot(:,end,end)];
                
                B_smoother = s5_B_scaling * Bc_ocean;
                B_smoother_inv = inv(B_smoother);
                %                                 if l_extra_ob  % Add in extra ob for smoother step
                %                                     for i=1:5
                %                                         if extra_ob_time(i) > 0
                %                                             ob_ix(i,extra_ob_time(i)) = 1;
                %                                         end
                %                                     end
                %                                 end

                for i_count_smoother = 1:s5_smoother_loops
                    if i_count_smoother == 1
                        u_lin = ua_plot_2;
                    else
                        u_lin = ua2_f;
                    end
                    innov_o = u_ob(:,1:nsteps*n_cycles_per_smoother) - u_lin(:,2:nsteps*n_cycles_per_smoother+1);
                    dX0_o=zeros(no,1);
                    [dXo_anal2, JXoInner2, dJXoInner2, ito2, rel_grad_o2] = minimize_mod_crit_NKN(dX0_o,'calcfg_ocean_l96c',max_iterations,tolerance,...
                        u_b,innov_o,u_lin,H,Boinv,Roinv,nsteps*n_cycles_per_smoother,h,na,no,Fx,Fy,alph,gamma,ob_ix);
                    %                         assim_scheme,h,nsteps*n_cycles_per_smoother,beta,rho,sigma, ...
                    %                         Omega,k,w_star,coupling_freq, u_b, B_smoother_inv, u_lin, 0, innov_o, var_ob5, ob_ix, l_fgat_s5);
                    if l_integration_coupled_s5 && i_smooth_iteration == s5_iterations && i_count_smoother == s5_smoother_loops
                        X_anal2 = u_lin(:,1) + [zeros(na,1);dXo_anal2];
                        X_temp = X_anal2;
                        for icycles = 1:n_cycles_per_smoother
                            if (update_method == 2 && icycles ~= 1)   % Adding the increment to the atmosphere forecast
                                X_temp(1:na) = X_temp(1:na) + Increment_Scaling*dXa_anal_cycle(:,icycles);
                            end
                            if (update_method == 2 && icycles == 1 )
                                X_temp(1:na) = u_lin(1:na,1);
                            end
                            za2_f = l96c_rk2(X_temp,h,nsteps,na,no,Fx,Fy,alph,gamma);
                            ua2_f(na+1:na+no,(icycles-1)*nsteps+1:icycles*nsteps) = za2_f(na+1:end,1:nsteps);
                            if (update_method == 2)
                                za2_f_new(:,(icycles-1)*nsteps+1:icycles*nsteps) = za2_f(:,1:nsteps);
                            end
                            if icycles == n_cycles_per_smoother
                                %                                         ua2_f(1:3,end) = [xa2_f(nsteps+1),ya2_f(nsteps+1),za2_f(nsteps+1)]';
                                ua2_f(na+1:na+no,end) = za2_f(na+1:na+no,nsteps+1);
                                %                                 ua2_f(1:na,end) = za2_f(1:na,nsteps+1);
                                ua2_f(1:na,:) = ua_plot_2(1:na,:);
                            end
                            %                             if (update_method == 1) %
                            X_temp(1:na) = ua_plot_2(1:na,icycles*nsteps+1); %
                            %                             elseif (update_method == 2 && reset_atmosphere == 0)
                            %                                 X_temp(1:na) = xa2_f(1:na,nsteps+1); %
                            %                             elseif (update_method == 2 && reset_atmosphere == 1)
                            %                                 X_temp(1:na) = ua_plot_2(1:na,icycles*nsteps+1); %
                            %                             end
                            X_temp(na+1:na+no) = za2_f(na+1:na+no,nsteps+1);
                        end
                    else
                        X_anal2 = u_lin(:,1) + [zeros(na,1);dXo_anal2];
                        [ya2_f] = l96c_rk2_ocean(ua_plot_2(1:na,1:end-1),X_anal2(na+1:end),...
                            h,nsteps*n_cycles_per_smoother,no,Fy,alph,gamma);
                        ua2_f(1:na,:)=ua_plot_2(1:na,:);
                        ua2_f(na+1:na+no,:)=ya2_f;
                    end
                end                
            end
            % ua2_f is the guess trajectory for the next iteration
            % (i_smooth_iteration), and is the final analysis if this
            % is the last iteration
            
            if l_lin_s5 == 0
                Bg_err_smoother((i_ob_pattern_repeats-1) * ...
                ob_pattern_repeat_freq + i_part_of_ob_pattern, :, :) = ua_plot_2(na+1:end,:) - u(na+1:end,:);
            end
            
            %                 if fcsteps > 0
            %                     [x_fc,y_fc,z_fc,w_fc,v_fc] = molteni_rk2(h,fcsteps,beta,rho,sigma,alpha, ...
            %                         Omega,k,w_star,ua2_f(1,end),ua2_f(2,end),ua2_f(3,end),ua2_f(4,end),ua2_f(5,end),coupling_freq);
            %                     ua2_fc = [x_fc y_fc z_fc w_fc v_fc]';
            %                     ua2_f = [ua2_f ua2_fc(:,2:end)];
            %                 end
            
        end
    end % i_smooth
    
    starttime = h * assim_steps * ((i_ob_pattern_repeats-1) * ob_pattern_repeat_freq + (i_part_of_ob_pattern-1));
    tvals=(starttime:h:starttime+h*assim_steps);

    %% plots of error norms
        if (l_plot_error_norm || l_plot_avg_error_norm)
            %% Store norm values for later plotting
            for i=1:n_cycles_per_smoother
                % Total norm
                bg_norm = vecnorm(rdivide(squeeze(ub_plot(:,i,:)) - u(:,(i-1)*nsteps+1:i*nsteps+1), SD'));
                anal_norm = vecnorm(rdivide(squeeze(ua_plot(:,i,:)) - u(:,(i-1)*nsteps+1:i*nsteps+1), SD'));
                Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 1) = bg_norm;
                Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 1) = anal_norm;
    
                % Atmospheric norm
                bg_norm = vecnorm(rdivide(squeeze(ub_plot(1:na,i,:)) - u(1:na,(i-1)*nsteps+1:i*nsteps+1), SD(1:na)'));
                anal_norm = vecnorm(rdivide(squeeze(ua_plot(1:na,i,:)) - u(1:na,(i-1)*nsteps+1:i*nsteps+1), SD(1:na)'));
                Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 2) = bg_norm;
                Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 2) = anal_norm;
    
                % Oceanic norm
                bg_norm = vecnorm(rdivide(squeeze(ub_plot(na+1:end,i,:)) - u(na+1:end,(i-1)*nsteps+1:i*nsteps+1), SD(na+1:end)'));
                anal_norm = vecnorm(rdivide(squeeze(ua_plot(na+1:end,i,:)) - u(na+1:end,(i-1)*nsteps+1:i*nsteps+1), SD(na+1:end)'));
                Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 3) = bg_norm;
                Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 3) = anal_norm;
            end
    
            if (assim_scheme == 5)
                % Total norm
                smoother_norm = vecnorm(rdivide(ua2_f - u, SD'));
                Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, :, 1) = smoother_norm;
    
                % Oceanic norm
                smoother_norm = vecnorm(rdivide(ua2_f(na+1:end,:) - u(na+1:end,:), SD(na+1:end)'));
                Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, :, 2) = smoother_norm;
            end
        end
        
        if l_plot_error_norm
            %% Plot evolution of error norms
            plot1a = figure(50);
            
            % Total norm
            subplot(3,1,1)
            hold on
            for i=1:n_cycles_per_smoother
                if i == 1
                    bg_norm_mean = 0;
                    anal_norm_mean = 0;
                end
                bg_norm_mean = bg_norm_mean + mean(0.5 * (Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 1:end-1, 1) + Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 2:end, 1)));
                anal_norm_mean = anal_norm_mean + mean(0.5 * (Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 1:end-1, 1) + Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 2:end, 1)));
                plot(tvals((i-1)*nsteps+1:i*nsteps+1), squeeze(Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 1)), 'Color', '#0072BD')
                plot(tvals((i-1)*nsteps+1:i*nsteps+1), squeeze(Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 1)), 'r')
                if i ~= n_cycles_per_smoother
                    xline(starttime+nsteps*h*i,'--','HandleVisibility','Off');
                elseif i_part_of_ob_pattern ~= ob_pattern_repeat_freq
                    xline(starttime+nsteps*h*i,'HandleVisibility','Off');
                else
                    xline(starttime+nsteps*h*i,'HandleVisibility','Off', 'Linewidth', 2);
                end
            end
            plot(tvals, bg_norm_mean / n_cycles_per_smoother * ones(length(tvals)), 'Color', '#0072BD', 'LineStyle', '--')
            plot(tvals, anal_norm_mean / n_cycles_per_smoother * ones(length(tvals)), '--r')
            if (assim_scheme == 5)
                plot(tvals, squeeze(Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, :, 1)), 'k')
                plot(tvals, mean(0.5 * (Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, 1:end-1, 1) + Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, 2:end, 1))) * ones(length(tvals)), '--k')
            end
            title('Error norms - Total')
            ylabel('Error norm')
            xlabel('Time')
        
            % Atmoshperic norm
            subplot(3,1,2)
            hold on
            for i=1:n_cycles_per_smoother
                if i == 1
                    bg_norm_mean = 0;
                    anal_norm_mean = 0;
                end
                bg_norm_mean = bg_norm_mean + mean(0.5 * (Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 1:end-1, 2) + Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 2:end, 2)));
                anal_norm_mean = anal_norm_mean + mean(0.5 * (Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 1:end-1, 2) + Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 2:end, 2)));
                plot(tvals((i-1)*nsteps+1:i*nsteps+1), squeeze(Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 2)), 'Color', '#0072BD')
                plot(tvals((i-1)*nsteps+1:i*nsteps+1), squeeze(Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 2)), 'r')
                if i ~= n_cycles_per_smoother
                    xline(starttime+nsteps*h*i,'--','HandleVisibility','Off');
                elseif i_part_of_ob_pattern ~= ob_pattern_repeat_freq
                    xline(starttime+nsteps*h*i,'HandleVisibility','Off');
                else
                    xline(starttime+nsteps*h*i,'HandleVisibility','Off', 'Linewidth', 2);
                end
            end
            plot(tvals, bg_norm_mean / n_cycles_per_smoother * ones(length(tvals)), 'Color', '#0072BD', 'LineStyle', '--')
            plot(tvals, anal_norm_mean / n_cycles_per_smoother * ones(length(tvals)), '--r')
            title('Error norms - Atmosphere')
            ylabel('Error norm')
            xlabel('Time')
        
            % Oceanic norm
            subplot(3,1,3)
            hold on
            for i=1:n_cycles_per_smoother
                if i == 1
                    bg_norm_mean = 0;
                    anal_norm_mean = 0;
                end
                bg_norm_mean = bg_norm_mean + mean(0.5 * (Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 1:end-1, 3) + Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 2:end, 3)));
                anal_norm_mean = anal_norm_mean + mean(0.5 * (Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 1:end-1, 3) + Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, 2:end, 3)));
                plot(tvals((i-1)*nsteps+1:i*nsteps+1), squeeze(Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 3)), 'Color', '#0072BD')
                plot(tvals((i-1)*nsteps+1:i*nsteps+1), squeeze(Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 3)), 'r')
                if i ~= n_cycles_per_smoother
                    xline(starttime+nsteps*h*i,'--','HandleVisibility','Off');
                elseif i_part_of_ob_pattern ~= ob_pattern_repeat_freq
                    xline(starttime+nsteps*h*i,'HandleVisibility','Off');
                else
                    xline(starttime+nsteps*h*i,'HandleVisibility','Off', 'Linewidth', 2);
                end
            end
            plot(tvals, bg_norm_mean / n_cycles_per_smoother * ones(length(tvals)), 'Color', '#0072BD', 'LineStyle', '--')
            plot(tvals, anal_norm_mean / n_cycles_per_smoother * ones(length(tvals)), '--r')
            if (assim_scheme == 5)
                plot(tvals, squeeze(Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, :, 2)), 'k')
                plot(tvals, mean(0.5 * (Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, 1:end-1, 2) + Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, 2:end, 2))) * ones(length(tvals)), '--k')
            end
            title('Error norms - Ocean')
            ylabel('Error norm')
            xlabel('Time')
        end
        l_marker = (i_ob_pattern_repeats-1)*assim_steps+(i_ob_pattern_repeats-1); 
        r_marker = l_marker+assim_steps;
        t_store(l_marker+1:r_marker) = (l_marker-1)*h:h:(r_marker-2)*h;
        ua2_f_store(:,l_marker+1:r_marker+1) = ua2_f(:,:);
        ua_plot_reshaped = reshape(permute(ua_plot(:,:,1:nsteps),[1,3,2]),[na+no,assim_steps]);
        ua_plot_reshaped_store(:,l_marker+1:r_marker) = ua_plot_reshaped(:,1:assim_steps);
        z_store(:,l_marker+1:r_marker+1) = z(:,:);
        u_ob_store(:,l_marker+1:r_marker) = u_ob(:,1:assim_steps);
end     % i_part_of_ob_pattern

u_ob_store(u_ob_store == 0) = nan;
ua2_f_store(ua2_f_store == 0) = nan;
ua_plot_reshaped_store(ua_plot_reshaped_store == 0) = nan;
index_show = 2;

plot_num2 = 20;
figure(plot_num2);
% analysis forecast of weakly coupled 4dvar:
plot(ua_plot_reshaped_store(index_show,:),'k-','DisplayName',strcat('Analysis Forecast of x',num2str(index_show))); hold on;
plot(z_store(index_show,:),'b-','DisplayName',strcat('True state of x',num2str(index_show))); hold on;
u_ob_plot = u_ob;
u_ob_plot(u_ob_plot == 0) = nan;
plot(u_ob_store(index_show,:),'r*','DisplayName',strcat('Observed state of x',num2str(index_show)))
xlabel('Assimilation Steps')
ylabel(strcat('x',num2str(index_show),'(t_{i})'))
legend show
plot_num2 = plot_num2 +1;

figure(plot_num2);
index_show = na+2;
% analysis forecast of weakly coupled 4dvar:
plot(ua_plot_reshaped_store(index_show,:),'k-','DisplayName',strcat('Analysis Forecast of y',num2str(index_show-na))); hold on;
plot(z_store(index_show,:),'b-','DisplayName',strcat('True State of y',num2str(index_show-na))); hold on;
plot(u_ob_store(index_show,:),'r*','DisplayName',strcat('Observed State of y',num2str(index_show-na)))
xlabel('Assimilation Steps')
ylabel(strcat('y',num2str(index_show-na),'(t_{i})'))
legend show
plot_num2 = plot_num2 +1;

figure(plot_num2);
index_show = 2+na;
% analysis forecast of weakly coupled 4dvar:
plot(ua_plot_reshaped_store(index_show,:),'g-','DisplayName',strcat('Analysis Forecast of y',num2str(index_show-na),',4dvar')); hold on;
plot(ua2_f_store(index_show,:),'k-','DisplayName',strcat('Model forecast of y',num2str(index_show-na),',after smoother')); hold on;...
plot(z_store(index_show,:),'b-','DisplayName',strcat('True state of y',num2str(index_show-na))); hold on;
plot(u_ob_store(index_show,:),'r*','DisplayName',strcat('Observed state of y',num2str(index_show-na)))
xlabel('Assimilation Steps')
ylabel(strcat('y',num2str(index_show-na),'(t_{i})'))
legend show
plot_num2 = plot_num2 + 1;

figure(plot_num2);
time_show = 200;
plot(ua_plot_reshaped_store(:,time_show),'k-','DisplayName',...
    strcat('Analysis Forecast of the State at Time Step t_{',num2str(time_show),'}'));
hold on;
plot(z_store(:,time_show),'b-','DisplayName',...
    strcat('True Value of the State at Time Step t_{',num2str(time_show),'}'))
hold on;
plot(u_ob_store(:,time_show),'r*','DisplayName',...
    strcat('Observed Value of the State at Time Step t_{',num2str(time_show),'}'))
xlabel('i')
ylabel('x_{i}')
legend show
plot_num2 = plot_num2 + 1;



