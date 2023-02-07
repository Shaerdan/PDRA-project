clc; clear all; close all;
tolerance = 1.0d-6;   % For solver
max_iterations = 20; % For solver
% Grad_Test = 1: turn on gradient tests for calcfg routines; = 0 turn off.
Grad_Test = 0;
%% rk2 solver parameters and model parameters for the l96 coupled model
nsteps = 50;
h=0.00125d0;
Fx=15;
Fy=8;
alph=0.5;
gamma= 0.6;
N = 40;
na = N; no = N; ntotal = na + no;
% loop controls:
outer_loops = 4;     % number of outerloop for weakly coupled standard 4dvar
s5_smoother_loops = 2;  % Number of outer loops for smoother step only
% n_cycles_per_smoother = 4; % number of short window cycles per smoother cycle
% n_ob_pattern_repeats = 4; % Total length of the run in terms of ob_pattern_repeat_freq
% method control:
min_method = 0; % 0 for NKN with Adjoint grad, 1 for fmincon with FD grad (bfgs)
min_method_smoother = 0; % smoother min method, same options as above
assim_scheme = 4;  % 5 for smoother method

if assim_scheme == 4
    n_cycles_per_smoother = 1;
    ob_pattern_repeat_freq = 4;
else
    n_cycles_per_smoother = 4;
    ob_pattern_repeat_freq = 1;
end

% this number of (smoother) assimilation cycles
assim_steps = nsteps*n_cycles_per_smoother;

l_fgat_s5 = 0;     % 0 = 4DVar for smoother step; 1 = 3DFGAT for smoother step
% data control:
l_newbg_xb = 1;
l_newobs = 1;
% plot control:
l_plot_convergence = 1;
l_plot_state = 0;
l_plot_error_norm = 1;
l_plot_avg_error_norm = 0;
l_plot_trajectories = 1;
%% Smoother setup
s5_B_scaling = 1;
s5_iterations = 1;
l_lin_s5 = 1;       % 0 = Take the analysis trajectory as both background and the first linearisation state;
l_integration_coupled_s5 = 1;   % At the last outer loop of the last iteration of the smoother step (which produces
% the final analysis), whether to integrate the smoother analysis at initial time using
% the coupled model; if yes, the atmospheric trajectory would be reset to that of the
% original atmospheric analyses at the beginning of each standard assimilation window
% before the coupled model integration is continued (this capability hasn't been designed
% for cases where IAU is used for the smoother step)
update_method = 1;
reset_atmosphere = 0;
Increment_Scaling = 1;
%% Setting up parameters for the assimilation

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

%% formulate background&observation error cov matrices, and H matrix:
number_of_samples = ntotal; % full sample size for (likely) nonsingular B
l_SpCov_SOAR = 1; % 0 for sampled covariance B, 1 for SOAR
L_atmos = 2; L_ocean = 4; % make these input variable
var_atmos_bg = 0.5; var_ocean_bg = 0.5;
[Bainv,Boinv,Ba,Bo,B,SD] = GetCovMatriceB(number_of_samples,h,assim_steps,na,no,Fx,Fy,alph,gamma,...
    l_SpCov_SOAR,L_atmos, L_ocean,var_atmos_bg, var_ocean_bg);
% B = blkdiag(Ba,Bo);
% observation pattern:
x_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
y_ob = (ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps:1:...
    ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
% observation stats:
var_ob = [1e-2, 1e-2];
R_atmos = var_ob(1)*eye(na,na); R_ocean = var_ob(2)*eye(no,no);
R = blkdiag(R_atmos,R_ocean);
Rinv = inv(R);
Rainv = inv(R_atmos);
Roinv = inv(R_ocean);
H = eye(ntotal,ntotal);

n_ob_pattern_repeats = 10;
for i_ob_pattern_repeats = 1:n_ob_pattern_repeats
    clear zb_f_chk za_chk zb_f_chk_store za_chk_store
    for i_part_of_ob_pattern = 1:ob_pattern_repeat_freq
        %% Start identical-twin set-up
        
        zb_plot=zeros(ntotal,n_cycles_per_smoother,nsteps+1);
        
        % Generate truth
        if (i_ob_pattern_repeats == 1 && i_part_of_ob_pattern == 1)
            x0_t = x0_init;
            y0_t = y0_init;
        else
            x0_t = x(:,assim_steps+1)';
            y0_t = y(:,assim_steps+1)';
        end
        z0_t = [x0_t,y0_t];
        [z] = l96c_rk2(z0_t,h,assim_steps,na,no,Fx,Fy,alph,gamma);
        x = z(1:na,:);
        y = z(na+1:end,:);
        z_t=[x0_t,y0_t]';
        X_t= z_t;
        
        % Generate background
        if (i_ob_pattern_repeats == 1 && i_part_of_ob_pattern == 1)
            if l_newbg_xb
                noise = randn(ntotal,1);
                z_b = z_t + sqrtm(B) * noise(1:ntotal);
                % debugging check
                figure(1000)
                plot((z_b-z_t)/mean(abs(z_t)))
                save(data_bgx_out,'z_b')
            else
                load(data_bgx_in,'z_b')
            end
        else
            z_b = za_plot(:,end,end);
        end
        % initialising za_plot
        za_plot=zeros(ntotal,n_cycles_per_smoother,nsteps+1);
        
        if l_newobs
            % Generate obs
            ob_ix=zeros(assim_steps,2);
            obs_noise=randn(ntotal,assim_steps);
            z_ob=zeros(ntotal,assim_steps);
            
            x_ob_local = x_ob((x_ob > (i_part_of_ob_pattern-1)*assim_steps) & (x_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
            y_ob_local = y_ob((y_ob > (i_part_of_ob_pattern-1)*assim_steps) & (y_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
            
            for i=x_ob_local
                ob_ix(i,1) = 1;
                z_ob(1:na,i) = H(1:na,1:na)*(z(1:na,i+1) + sqrt(var_ob(1))*obs_noise(1:na,i));
            end
            for i=y_ob_local
                ob_ix(i,2) = 1;
                z_ob(na+1:ntotal,i) =  H(na+1:ntotal,na+1:ntotal)*(z(na+1:ntotal,i+1) + sqrt(var_ob(2))*obs_noise(na+1:ntotal,i));
            end
            save(data_obs_out,'z_ob','ob_ix')
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
                    X_b=z_b;
                    X_t=z_t;
                    Errormat_Tot_Atm_Oce_anal= zeros(s5_iterations, outer_loops, 3);
                    % h,nsteps+fcsteps,beta,rho,sigma,alpha_b, ...
                    %    Omega,k,w_star,z_b(1),z_b(2),z_b(3),z_b(4),z_b(5),coupling_freq
                    %                     [zb_f,M] = transmat_lor96c(z_b,nsteps*h,h,na,Fx,Fy,alph,gamma);
                    [zb_f] = l96c_rk2(z_b,h,nsteps,na,no,Fx,Fy,alph,gamma);
                    xb_f = zb_f(1:na,:);
                    yb_f = zb_f(na+1:ntotal,:);
                    zb_plot(:,i_cycles,:) = zb_f;
                    zb_f_chk(:,(i_cycles-1)*nsteps+1:i_cycles*nsteps) = zb_f(:,1:end-1);
                end
                
                %% Perform assimilation of the weakly coupled standard 4DVAR:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i_count = 1:outer_loops
                    % Calculate innovations using coupled model
                    if (i_smooth_iteration == 1 && i_count == 1)
                        % Use background trajectory
                        z0 = z_b;
                        x_lin = xb_f(:,1:nsteps+1);
                        y_lin = yb_f(:,1:nsteps+1);
                    else
                        % Calculate nonlinear trajectory with coupled model
                        [zb_f] = l96c_rk2(z0,h,nsteps,na,no,Fx,Fy,alph,gamma);
                    end
                    z_lin = zb_f;
                    %                     z_lin = [x_lin,y_lin,z_lin,w_lin,v_lin]';
                    innov = z_ob(:,(i_cycles-1)*nsteps+1:i_cycles*nsteps) - H*z_lin(:,2:nsteps+1); % Assume no ob at time zero.
                    % Includes where there are no obs, but can use ob_ix
                    % Perform separate minimisations
                    
                    %% weakly coupled standard 4dvar, minimisation routine:
                    if min_method == 1
                        % Atmosphere
                        options0 = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',false,...
                            'PlotFcn','optimplotfval','MaxIterations',20);
                        dX0_a=zeros(na,1);
                        [dXa_anal,JXaInner,exitflag1,output1,lambda1,dJXaInner,hessian1] = fmincon(@(Xmin) calcfg_atmos_l96c(Xmin,zb_plot(:,i_cycles,1),innov,z_lin,H,...
                            Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix(:,1)),dX0_a,[],[],[],[],[],[],[],options0);
                        ita = output1.iterations;
                        % Ocean
                        dX0_o=zeros(no,1);
                        [dXo_anal,JXoInner,exitflag2,output2,lambda2,dJXoInner,hessian2] = fmincon(@(Xmin) calcfg_ocean_l96c(Xmin,zb_plot(:,i_cycles,1),innov,z_lin,H,...
                            Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix(:,2)),dX0_o);
                        ito = output2.iterations;
                        
                    elseif min_method == 0
                        % Atmosphere
                        dX0_a=zeros(na,1);
                        [dXa_anal, JXaInner, dJXaInner, ita, rel_grad_a] = minimize_mod_crit_NKN(dX0_a,'calcfg_atmos_l96c',max_iterations,tolerance, ...
                            zb_plot(:,i_cycles,1),innov,z_lin,H,...
                            Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix(:,1));
                        
                        % Ocean
                        dX0_o=zeros(no,1);
                        [dXo_anal, JXoInner, dJXoInner, ito, rel_grad_o] = minimize_mod_crit_NKN(dX0_o,'calcfg_ocean_l96c',max_iterations,tolerance, ...
                            zb_plot(:,i_cycles,1),innov,z_lin,H,...
                            Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix(:,2));
                        format short
                        IT_Outer_ITa_Ito = [i_count, ita, ito];
                        format longe
                        rel_grad_a;
                        rel_grad_o;
                        Cost_RG_Atmos = [JXaInner, rel_grad_a];
                        Cost_RG_Ocean = [JXoInner, rel_grad_o];
                        if Grad_Test == 1
                            gradient_test(rand(na,1),rand(no,1),zb_plot(:,i_cycles,1),innov,...
                                z_lin,H,ob_ix,Bainv,Rainv,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma)
                        end
                    end % end of minimisation routines of weakly coupled
                    %%%%%%%%%%%%
                    
                    % Update fields
                    delta_z0 = [dXa_anal; dXo_anal];
                    z0 = z0 + delta_z0;
                    
                    %%%%%%%%%%%
                    format long
                    X_anal = z0;
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
                    
                end % end outerloop
                %                     dXa_anal_cycle(:,i_cycles) = z0(1:na) - X_b(1:na);
                %                     Count_dXa = Count_dXa + 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %% Display error norms and store background errors for diagnostics
                disp('***Final Results***')
                format long
                Xb_Xt_Xanal =[X_b(1), X_t(1), X_anal(1)]
                
                %%%%%%%%%%%%
                if min_method == 0
                    plot_convergence(JXa,dJXa,JXo,dJXo,i_cycles,i_part_of_ob_pattern,n_cycles_per_smoother,...
                        ob_pattern_repeat_freq,l_plot_convergence,assim_scheme)
                end
                
                %% Run forecast
                xa=X_anal(1:na);
                ya=X_anal(na+1:ntotal);
                
                za = [xa;ya];
                [za_f] = l96c_rk2(za,h,nsteps,na,no,Fx,Fy,alph,gamma);
                
                za_plot(:,i_cycles,:) = za_f;
                za_chk(:,(i_cycles-1)*nsteps+1:(i_cycles-1)*nsteps+nsteps) = za_f(:,1:end-1);
                z_ob_chk = z_ob;
                z_ob_chk(z_ob_chk == 0) = nan;
                if assim_scheme == 4
                    indx_str = (i_part_of_ob_pattern-1)*nsteps+1;
                    indx_end = (i_part_of_ob_pattern-1)*nsteps + nsteps;
                    za_chk_store(:,indx_str:indx_end) = za_chk;
                    z_store(:,indx_str:indx_end+1) = z;
                    zb_f_chk_store(:,indx_str:indx_end) = zb_f_chk;
                    z_ob_chk_store(:,indx_str:indx_end) = z_ob_chk;
                    error_norm_analysis(indx_str:indx_end) = vecnorm(za_chk - z(:,1:end-1));
                    error_norm_bg(indx_str:indx_end) = vecnorm(zb_f_chk - z(:,1:end-1));
                else
                    indx_str = (i_cycles-1)*nsteps+1;
                    indx_end = (i_cycles-1)*nsteps + nsteps;
                    error_norm_analysis(indx_str:indx_end) = vecnorm(za_chk(2,indx_left:indx_right)...
                                                            -z(2,1:end-1));
                    error_norm_bg(indx_str:indx_end) = vecnorm(zb_f_chk(2,indx_left:indx_right)...
                        - z(2,1:end-1));
                end
                % debug plotting za_f vs z;
                indx_show = 2;
                if i_cycles == n_cycles_per_smoother && assim_scheme == 5
                    figure(200 + i_ob_pattern_repeats)
                    plot(za_chk(indx_show,:),'k-*','DisplayName','Analysis Forecast'); hold on;
                    plot(z(indx_show,:),'r-','DisplayName','Ground Truth'); hold on;
                    plot(zb_f_chk(indx_show,:),'b-','DisplayName','Background Forecast'); hold on;
                    plot(z_ob_chk(indx_show,:),'go','DisplayName','Observation');
                    xlabel('Assimilation Steps')
                    legend show
                    figure(1200 + i_ob_pattern_repeats)
                    plot(error_norm_bg,'b-*','DisplayName','Background Trajectory Error Norm'); hold on;
                    plot(error_norm_analysis,'k-*','DisplayName','Analysis Trajectory Error Norm')
                    legend show
                elseif i_part_of_ob_pattern == ob_pattern_repeat_freq && assim_scheme == 4
                    figure(200 + i_ob_pattern_repeats)
                    plot(za_chk_store(indx_show,:),'k-*','DisplayName','Analysis Forecast'); hold on;
                    plot(z_store(indx_show,:),'r-','DisplayName','Ground Truth'); hold on;
                    plot(zb_f_chk_store(indx_show,:),'b-','DisplayName','Background Forecast'); hold on;
                    plot(z_ob_chk_store(indx_show,:),'go','DisplayName','Observation');
                    xlabel('Assimilation Steps')
                    legend show
                    figure(1200 + i_ob_pattern_repeats)
                    plot(error_norm_bg,'b-*','DisplayName','Background Trajectory Error Norm'); hold on;
                    plot(error_norm_analysis,'k-*','DisplayName','Analysis Trajectory Error Norm')
                    legend show                  
                end
                if i_smooth_iteration == 1
                    %% Set background and truth for next cycle
                    z_b = za_f(:,nsteps+1);
                    z_t = z(:,i_cycles*nsteps+1);
                end
            end % i_cycles
            if assim_scheme == 5
                [za2_f,z0] = smoother_step(za_plot,zb_plot,assim_steps,s5_B_scaling,...
                    Bo,Roinv,H,s5_smoother_loops,z_ob,n_cycles_per_smoother,l_integration_coupled_s5,...
                    s5_iterations,i_smooth_iteration,...
                    h,nsteps,na,no,Fx,Fy,alph,gamma,ob_ix,i_ob_pattern_repeats,ob_pattern_repeat_freq,...
                    i_part_of_ob_pattern,l_lin_s5,max_iterations,tolerance,min_method_smoother);
                figure(400 + i_ob_pattern_repeats)
                plot(za2_f(2+na,:),'b-o','DisplayName','PostSmoother Analysis Forecast'); hold on;
                plot(za_chk(2+na,:),'k-*','DisplayName','Presmoother Analysis Forecast'); hold on;
                plot(z(2+na,:),'r-','DisplayName','True State');
                xlabel('Assimilation Steps')
                legend show
            end
            
        end
    end % i_smooth
    
    %% insert plotting module here
    %     starttime = h * assim_steps * ((i_ob_pattern_repeats-1) * ob_pattern_repeat_freq + (i_part_of_ob_pattern-1));
    %     tvals=(starttime:h:starttime+h*assim_steps);
    
    %%
    %     ua_plot_reshaped = reshape(permute(za_plot(:,:,1:nsteps),[1,3,2]),[na+no,assim_steps]);
    %     ua_plot_reshaped_store(:,l_marker+1:r_marker) = ua_plot_reshaped(:,1:assim_steps);
    %     z_store(:,l_marker+1:r_marker+1) = z(:,:);
    %     u_ob_store(:,l_marker+1:r_marker) = z_ob(:,1:assim_steps);
end     % i_part_of_ob_pattern
save_all_figures = 1;
if save_all_figures == 1
 figHandles = findall(0,'Type','figure');
 for i = 1:numel(figHandles)
     fn = tempname(strcat('C:\06022023\results\'));  %in this example, we'll save to a temp directory.
     export_fig(fn, '-png', figHandles(i))
 end
end
% u_ob_store(u_ob_store == 0) = nan;
% if assim_scheme == 5
%     ua2_f_store(ua2_f_store == 0) = nan;
% end
% ua_plot_reshaped_store(ua_plot_reshaped_store == 0) = nan;
% index_show = 2;
%
% plot_num2 = 20;
% figure(20);
% % analysis forecast of weakly coupled 4dvar:
% plot(ua_plot_reshaped_store(index_show,:),'k-','DisplayName',strcat('Analysis Forecast of x',num2str(index_show))); hold on;
% plot(z_store(index_show,:),'b-','DisplayName',strcat('True state of x',num2str(index_show))); hold on;
% u_ob_plot = z_ob;
% u_ob_plot(u_ob_plot == 0) = nan;
% plot(u_ob_store(index_show,:),'r*','DisplayName',strcat('Observed state of x',num2str(index_show)))
% xlabel('Assimilation Steps')
% ylabel(strcat('x',num2str(index_show),'(t_{i})'))
% legend show
% plot_num2 = plot_num2 +1;
%
% figure(21);
% index_show = na+2;
% % analysis forecast of weakly coupled 4dvar:
% plot(ua_plot_reshaped_store(index_show,:),'k-','DisplayName',strcat('Analysis Forecast of y',num2str(index_show-na))); hold on;
% plot(z_store(index_show,:),'b-','DisplayName',strcat('True State of y',num2str(index_show-na))); hold on;
% plot(u_ob_store(index_show,:),'r*','DisplayName',strcat('Observed State of y',num2str(index_show-na)))
% xlabel('Assimilation Steps')
% ylabel(strcat('y',num2str(index_show-na),'(t_{i})'))
% legend show
% plot_num2 = plot_num2 +1;
%
% if assim_scheme == 5
%     figure(22);
%     index_show = 2+na;
%     % analysis forecast of weakly coupled 4dvar:
%     plot(ua_plot_reshaped_store(index_show,:),'g-','DisplayName',strcat('Analysis Forecast of y',num2str(index_show-na),',4dvar')); hold on;
%     plot(ua2_f_store(index_show,:),'k-','DisplayName',strcat('Model forecast of y',num2str(index_show-na),',after smoother')); hold on;...
%         plot(z_store(index_show,:),'b-','DisplayName',strcat('True state of y',num2str(index_show-na))); hold on;
%     plot(u_ob_store(index_show,:),'r*','DisplayName',strcat('Observed state of y',num2str(index_show-na)))
%     xlabel('Assimilation Steps')
%     ylabel(strcat('y',num2str(index_show-na),'(t_{i})'))
%     legend show
%     plot_num2 = plot_num2 + 1;
% end
%
% figure(23);
% time_show = 200;
% plot(ua_plot_reshaped_store(:,time_show),'k-','DisplayName',...
%     strcat('Analysis Forecast of the State at Time Step t_{',num2str(time_show),'}'));
% hold on;
% plot(z_store(:,time_show),'b-','DisplayName',...
%     strcat('True Value of the State at Time Step t_{',num2str(time_show),'}'))
% hold on;
% plot(u_ob_store(:,time_show),'r*','DisplayName',...
%     strcat('Observed Value of the State at Time Step t_{',num2str(time_show),'}'))
% xlabel('i')
% ylabel('x_{i}')
% legend show
% plot_num2 = plot_num2 + 1;

 


