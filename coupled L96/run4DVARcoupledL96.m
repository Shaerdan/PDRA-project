clc;
clear all;
close all;

n_repeat_no_cycle = 1;
%% rk2 solver parameters and model parameters for the l96 coupled model
for i_trial = 1:2
    rng(1)  % Fix random seed
    tolerance = 1.0d-6;   % For solver
    max_iterations = 100; % For solver
    % Grad_Test = 1: turn on gradient tests for calcfg routines; = 0 turn off.
    Grad_Test = 0;
    save_all_figures = 0;
    dirname = strcat('C:\results\DA\14022023\8randonecycle\');
    nsteps = 4;
    h=0.0125d0;
    Fx=15;
    Fy=8;
    alph=0.5;
    gamma= 0;
    N = 40;
    na = N; no = N; ntotal = na + no;
    var_atmos_bg = 1e-0; var_ocean_bg = 1e-0;
    var_ob = [1e-0, 1e-0];
    % loop controls:
    n_ob_pattern_repeats = 5000;
    outer_loops = 1;     % number of outerloop for weakly coupled standard 4dvar
    s5_smoother_loops = 2;  % Number of outer loops for smoother step only
    % method control:
    min_method = 0; % 0 for NKN with Adjoint grad, 1 for fmincon with FD grad (bfgs)
    min_method_smoother = 0; % smoother min method, same options as above
    schemes_trial = [4 5];
    assim_scheme = schemes_trial(i_trial);  % 5 for smoother method
    
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
    l_plot_error_norm = 0;
    l_plot_avg_error_norm = 1;
    l_plot_avg_error_norm_compare = 1;
    l_plot_trajectories = 0;
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
    update_method = 2;
    Increment_Scaling = 1;
    %% Setting up parameters for the assimilation
    
    data_bgx_out='data_bgBx.mat';  % Output file for xb
    data_obs_out='data_obs.mat';   % Output file or x_obs
    
    xvals=1:na; % atmosphere grid indeces
    yvals=1:no; % ocean grid indeces
    for i_repeat_one_cycle = 1: n_repeat_no_cycle
        
        x0_init=sin(xvals/(na-1)*2*pi);
        y0_init=cos(5*yvals/(no-1)*2*pi);
        % fun experiments with the model for N<=4:
        if n_repeat_no_cycle == 1
            i_model_tn = randperm(1000,10) + 200;
            [z_chk] = l96c_rk2([x0_init';y0_init'],h,i_model_tn(end)*assim_steps,na,no,Fx,Fy,alph,gamma);
            x0_init = z_chk(1:na,end)';
            y0_init = z_chk(na+1:end,end)';
%             figure(900)
%             plot3(z_chk(1+na,:),z_chk(2+na,:),z_chk(3+na,:),'r-'); hold on; ...
%                 plot3(z_chk(1,:),z_chk(2,:),z_chk(3,:),'k-');
        end
        %% formulate background&observation error cov matrices, and H matrix:
        number_of_samples = ntotal; % full sample size for (likely) nonsingular B
        l_SpCov_SOAR = 1; % 0 for sampled covariance B, 1 for SOAR
        L_atmos = 4; L_ocean = 2; % make these input variable
        
        [Bainv,Boinv,Ba,Bo,B,SD] = GetCovMatriceB(number_of_samples,h,assim_steps,na,no,Fx,Fy,alph,gamma,...
            l_SpCov_SOAR,L_atmos, L_ocean,var_atmos_bg, var_ocean_bg);
%         Bainv = (var_atmos_bg^(-1))*eye(na,na); Boinv = (var_atmos_bg^(-1))*eye(no,no);
        % B = blkdiag(Ba,Bo);
        % observation pattern:
        x_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
        y_ob = (ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps:1:...
            ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
        % observation stats:
        R_atmos = var_ob(1)*eye(na,na); R_ocean = var_ob(2)*eye(no,no);
        R = blkdiag(R_atmos,R_ocean);
        Rinv = inv(R);
        Rainv = inv(R_atmos);
        Roinv = inv(R_ocean);
        space_skip = 1;
        H_space_pattern = 1:space_skip:ntotal;
        H_diag = zeros(1,ntotal);
        H_diag(H_space_pattern) = 1;
        H = diag(H_diag);
        N_Obs_Num_Spatial = length(H_space_pattern);
        N_Obs_Num_Spatial_o = N_Obs_Num_Spatial/2;
        N_Obs_Num_Spatial_a = N_Obs_Num_Spatial_o;
        % H = eye(ntotal,ntotal);
        %     basetime = 0;
        
        %   initialise error norm matrices:
        if l_plot_avg_error_norm_compare
            if n_repeat_no_cycle ~=1
                n_store_norm = n_repeat_no_cycle;
            else
                n_store_norm = n_ob_pattern_repeats;
            end
            if i_trial == 1
                Err_norm_bg_1 = zeros(n_store_norm, ob_pattern_repeat_freq, n_cycles_per_smoother, nsteps + 1, 3);      % Last dimension: Total, Atmosphere, Ocean
                Err_norm_anal_1 = zeros(n_store_norm, ob_pattern_repeat_freq, n_cycles_per_smoother, nsteps + 1, 3);    % Last dimension: Total, Atmosphere, Ocean
            elseif i_trial == 2
                Err_norm_bg_2 = zeros(n_store_norm, ob_pattern_repeat_freq, n_cycles_per_smoother, nsteps + 1, 3);      % Last dimension: Total, Atmosphere, Ocean
                Err_norm_anal_2 = zeros(n_store_norm, ob_pattern_repeat_freq, n_cycles_per_smoother, nsteps + 1, 3);    % Last dimension: Total, Atmosphere, Ocean
                Err_norm_smoother = zeros(n_store_norm, ob_pattern_repeat_freq, nsteps * n_cycles_per_smoother + 1, 2);    % Last dimension: Total, Ocean
            end
        end
        
        Bg_err = zeros(n_ob_pattern_repeats * ob_pattern_repeat_freq * n_cycles_per_smoother, 5);
        
        for i_ob_pattern_repeats = 1:n_ob_pattern_repeats
            %         clear zb_f_chk za_chk zb_f_chk_store za_chk_store
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
                % test running the model sequentially
                %         z_test1 = l96c_rk2(z0_t,h,nsteps,na,no,Fx,Fy,alph,gamma);
                %         z_test2 = l96c_rk2(z_test1(:,end),h,nsteps,na,no,Fx,Fy,alph,gamma);
                %         z_test3 = l96c_rk2(z_test2(:,end),h,nsteps,na,no,Fx,Fy,alph,gamma);
                %         z_test4 = l96c_rk2(z_test3(:,end),h,nsteps,na,no,Fx,Fy,alph,gamma);
                
                x = z(1:na,:);
                y = z(na+1:end,:);
                z_t=[x0_t,y0_t]'; % truth at t_{(k-1)*assim_steps+1}, k is the repeated obs pattern
                X_t= z_t;
                
                % Generate background
                if (i_ob_pattern_repeats == 1 && i_part_of_ob_pattern == 1)
                    if l_newbg_xb
                        noise = randn(ntotal,1);
                        z_b = z_t + sqrtm(B) * noise(1:ntotal);
                        save(data_bgx_out,'z_b')
                    else
                        %load(data_bgx_in,'z_b')
                    end
                elseif assim_scheme == 5
                    z_b = za2_f(:,assim_steps+1);
                elseif assim_scheme == 4
                    z_b = za_plot(:,end,end);
                end
                disp(strcat('Pattern Repeats = ',num2str(i_ob_pattern_repeats)))
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
                            
                            [zb_f] = l96c_rk2(z_b,h,nsteps,na,no,Fx,Fy,alph,gamma);
                            xb_f = zb_f(1:na,:);
                            yb_f = zb_f(na+1:ntotal,:);
                            zb_plot(:,i_cycles,:) = zb_f;
%                             zb_f_chk(:,(i_cycles-1)*nsteps+1:i_cycles*nsteps+1) = zb_f;
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
                                z_lin = [x_lin;y_lin];
                            else
                                % Calculate nonlinear trajectory with coupled model
                                [z_lin] = l96c_rk2(z0,h,nsteps,na,no,Fx,Fy,alph,gamma);
                            end
                            %                     z_lin = [x_lin,y_lin,z_lin,w_lin,v_lin]';
                            innov = z_ob(:,(i_cycles-1)*nsteps+1:i_cycles*nsteps) - H*z_lin(:,2:nsteps+1); % Assume no ob at time zero.
                            % Includes where there are no obs, but can use ob_ix
                            % Perform separate minimisations
                            
                            %% weakly coupled standard 4dvar, minimisation routine:
                            if min_method == 1
                                % Atmosphere
                                options0 = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',false,...
                                    'PlotFcn','optimplotfval','MaxIterations',max_iterations);
                                dX0_a=zeros(na,1);
                                [dXa_anal,JXaInner,exitflag1,output1,lambda1,dJXaInner,hessian1] = fmincon(@(Xmin) calcfg_atmos_l96c(Xmin,zb_plot(:,i_cycles,1),innov,z_lin,H,...
                                    Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix(:,1)),dX0_a,[],[],[],[],[],[],[],options0);
                                ita = output1.iterations;
                                % Ocean
                                dX0_o=zeros(no,1);
                                [dXo_anal,JXoInner,exitflag2,output2,lambda2,dJXoInner,hessian2] = fmincon(@(Xmin) calcfg_ocean_l96c(Xmin,zb_plot(:,i_cycles,1),innov,z_lin,H,...
                                    Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix(:,2)),dX0_o);
                                ito = output2.iterations;
                                
                            elseif min_method == 0
                                % Atmosphere
                                dX0_a=zeros(na,1);
                                [dXa_anal, JXaInner, dJXaInner, ita, rel_grad_a] = minimize_mod_crit_NKN(dX0_a,'calcfg_atmos_l96c',max_iterations,tolerance, ...
                                    zb_plot(:,i_cycles,1),innov,z_lin,H,...
                                    Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix(:,1));
                                
                                % Ocean
                                dX0_o=zeros(no,1);
                                [dXo_anal, JXoInner, dJXoInner, ito, rel_grad_o] = minimize_mod_crit_NKN(dX0_o,'calcfg_ocean_l96c',max_iterations,tolerance, ...
                                    zb_plot(:,i_cycles,1),innov,z_lin,H,...
                                    Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix(:,2));
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
                        
                        if i_smooth_iteration == 1
                            %% Set background and truth for next cycle
                            z_b = za_f(:,nsteps+1);
                            z_t = z(:,i_cycles*nsteps+1);
                        end
                    end % i_cycles
                    
                    if assim_scheme == 5
                        [za2_f,z0] = smoother_step(za_plot,zb_plot,assim_steps,s5_B_scaling,...
                            Bo,Roinv,H,s5_smoother_loops,z_ob,n_cycles_per_smoother,...
                            h,nsteps,na,no,Fx,Fy,alph,gamma,ob_ix,l_lin_s5,...
                            max_iterations,tolerance,min_method_smoother);
                    end
                end % i_smooth
                %%%%%%%%% plotting module storage:
                if l_plot_avg_error_norm_compare
                    %% Store norm values for later plotting
                    if n_repeat_no_cycle ~=1
                        i_store_norm = i_repeat_one_cycle;
                    else
                        i_store_norm = i_ob_pattern_repeats;
                    end
                    for i=1:n_cycles_per_smoother
                        % Total norm
                        bg_norm = vecnorm(rdivide(squeeze(zb_plot(:,i,:)) - z(:,(i-1)*nsteps+1:i*nsteps+1), N_Obs_Num_Spatial));
                        anal_norm = vecnorm(rdivide(squeeze(za_plot(:,i,:)) - z(:,(i-1)*nsteps+1:i*nsteps+1), N_Obs_Num_Spatial));
                        if i_trial == 1
                            Err_norm_bg_1(i_store_norm, i_part_of_ob_pattern, i, :, 1) = bg_norm;
                            Err_norm_anal_1(i_store_norm, i_part_of_ob_pattern, i, :, 1) = anal_norm;
                        elseif i_trial == 2
                            Err_norm_bg_2(i_store_norm, i_part_of_ob_pattern, i, :, 1) = bg_norm;
                            Err_norm_anal_2(i_store_norm, i_part_of_ob_pattern, i, :, 1) = anal_norm;
                        end
                        
                        % Atmospheric norm
                        bg_norm = vecnorm(rdivide(squeeze(zb_plot(1:na,i,:)) - z(1:na,(i-1)*nsteps+1:i*nsteps+1), N_Obs_Num_Spatial_a));
                        anal_norm = vecnorm(rdivide(squeeze(za_plot(1:na,i,:)) - z(1:na,(i-1)*nsteps+1:i*nsteps+1), N_Obs_Num_Spatial_a));
                        if i_trial == 1
                            Err_norm_bg_1(i_store_norm, i_part_of_ob_pattern, i, :, 2) = bg_norm;
                            Err_norm_anal_1(i_store_norm, i_part_of_ob_pattern, i, :, 2) = anal_norm;
                        elseif i_trial == 2
                            Err_norm_bg_2(i_store_norm, i_part_of_ob_pattern, i, :, 2) = bg_norm;
                            Err_norm_anal_2(i_store_norm, i_part_of_ob_pattern, i, :, 2) = anal_norm;
                        end
                        
                        % Oceanic norm
                        bg_norm = vecnorm(rdivide(squeeze(zb_plot(na+1:ntotal,i,:)) - z(na+1:ntotal,(i-1)*nsteps+1:i*nsteps+1), N_Obs_Num_Spatial_o));
                        anal_norm = vecnorm(rdivide(squeeze(za_plot(na+1:ntotal,i,:)) - z(na+1:ntotal,(i-1)*nsteps+1:i*nsteps+1), N_Obs_Num_Spatial_o));
                        if i_trial == 1
                            Err_norm_bg_1(i_store_norm, i_part_of_ob_pattern, i, :, 3) = bg_norm;
                            Err_norm_anal_1(i_store_norm, i_part_of_ob_pattern, i, :, 3) = anal_norm;
                        elseif i_trial == 2
                            Err_norm_bg_2(i_store_norm, i_part_of_ob_pattern, i, :, 3) = bg_norm;
                            Err_norm_anal_2(i_store_norm, i_part_of_ob_pattern, i, :, 3) = anal_norm;
                        end
                    end
                    
                    if assim_scheme == 5 && i_trial == 2
                        % Total norm
                        smoother_norm = vecnorm(rdivide(za2_f - z, N_Obs_Num_Spatial));
                        Err_norm_smoother(i_store_norm, i_part_of_ob_pattern, :, 1) = smoother_norm;
                        
                        % Oceanic norm
                        smoother_norm = vecnorm(rdivide(za2_f(na+1:ntotal,:) - z(na+1:ntotal,:), N_Obs_Num_Spatial_o));
                        Err_norm_smoother(i_store_norm, i_part_of_ob_pattern, :, 2) = smoother_norm;
                    end
                end % plotting module storage end
                
            end  % i_part_of_ob_pattern
            disp('long window end')
        end     % i_ob_pattern_repeats
    end % one cycle repeat trials (no long cycling)
    
    %%%%%%%%%%% Plotting module
    if l_plot_avg_error_norm && i_repeat_one_cycle == n_repeat_no_cycle
        %% Plot error norms averaged over (smoother) cycles
        if i_trial == 1
            plot1a = figure(50);
        elseif i_trial == 2
            plot1b = figure(55);
        end
        
        % Total norm
        subplot(3,1,1)
        hold on
        if i_trial == 1
            for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
                basetime = (i-1) * nsteps * h;
                bg_error = squeeze(mean(Err_norm_bg_1(:,i,:,:,1), 1));
                plot((basetime:h:basetime+nsteps*h), bg_error, 'Color', '#0072BD')
                anal_error = squeeze(mean(Err_norm_anal_1(:,i,:,:,1), 1));
                plot((basetime:h:basetime+nsteps*h), anal_error, 'r')
                xline(basetime + nsteps * h,'HandleVisibility','Off');
            end
        elseif i_trial == 2
            for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
                basetime = (i-1) * nsteps * h;
                bg_error = squeeze(mean(Err_norm_bg_2(:,:,i,:,1), 1));
                plot((basetime:h:basetime+nsteps*h), bg_error, 'Color', '#0072BD')
                xline(basetime + nsteps * h,'HandleVisibility','Off');
            end
            anal_error = squeeze(mean(Err_norm_smoother(:,:,:,1), 1));
            plot((0:h:nsteps*n_cycles_per_smoother*h), anal_error, 'r')
        end
        title(sprintf('Error norms averaged over %i windows of an identical observation pattern - Total', i_store_norm))
        ylabel('Error norm')
        xlabel('Time within the window')
        hold off
        
        % Atmospheric norm
        subplot(3,1,2)
        hold on
        if i_trial == 1
            for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
                basetime = (i-1) * nsteps * h;
                bg_error = squeeze(mean(Err_norm_bg_1(:,i,:,:,2), 1));
                plot((basetime:h:basetime+nsteps*h), bg_error, 'Color', '#0072BD')
                anal_error = squeeze(mean(Err_norm_anal_1(:,i,:,:,2), 1));
                plot((basetime:h:basetime+nsteps*h), anal_error, 'r')
                xline(basetime + nsteps * h,'HandleVisibility','Off');
            end
        elseif i_trial == 2
            for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
                basetime = (i-1) * nsteps * h;
                bg_error = squeeze(mean(Err_norm_bg_2(:,:,i,:,2), 1));
                plot((basetime:h:basetime+nsteps*h), bg_error, 'Color', '#0072BD')
                anal_error = squeeze(mean(Err_norm_anal_2(:,:,i,:,2), 1));
                plot((basetime:h:basetime+nsteps*h), anal_error, 'r')
                xline(basetime + nsteps * h,'HandleVisibility','Off');
            end
        end
        title(sprintf('Error norms averaged over %i windows of an identical observation pattern - Atmoshpere', i_store_norm))
        ylabel('Error norm')
        xlabel('Time within the window')
        hold off
        
        % Ocenaic norm
        subplot(3,1,3)
        hold on
        if i_trial == 1
            for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
                basetime = (i-1) * nsteps * h;
                bg_error = squeeze(mean(Err_norm_bg_1(:,i,:,:,3), 1));
                plot((basetime:h:basetime+nsteps*h), bg_error, 'Color', '#0072BD')
                anal_error = squeeze(mean(Err_norm_anal_1(:,i,:,:,3), 1));
                plot((basetime:h:basetime+nsteps*h), anal_error, 'r')
                xline(basetime + nsteps * h,'HandleVisibility','Off');
            end
        elseif i_trial == 2
            for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
                basetime = (i-1) * nsteps * h;
                bg_error = squeeze(mean(Err_norm_bg_2(:,:,i,:,3), 1));
                plot((basetime:h:basetime+nsteps*h), bg_error, 'Color', '#0072BD')
                xline(basetime + nsteps * h,'HandleVisibility','Off');
            end
            anal_error = squeeze(mean(Err_norm_smoother(:,:,:,2), 1));
            plot((0:h:nsteps*n_cycles_per_smoother*h), anal_error, 'r')
        end
        title(sprintf('Error norms averaged over %i windows of an identical observation pattern - Ocean', i_store_norm))
        ylabel('Error norm')
        xlabel('Time within the window')
        hold off
    end  % end of plotting module
end % i_trial
if l_plot_avg_error_norm_compare && i_trial == 2
    %% Plot error norms averaged over (smoother) cycles
    plot2 = figure(60);
    
    % Total norm
    subplot(3,1,1)
    hold on
    for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
        basetime = (i-1) * nsteps * h;
        bg_error_compare = (squeeze(mean(Err_norm_bg_2(:,:,i,:,1), 1)) ...
            - squeeze(mean(Err_norm_bg_1(:,i,:,:,1), 1))) ...
            ./ squeeze(mean(Err_norm_bg_1(:,i,:,:,1), 1));
        if i ~= n_cycles_per_smoother*ob_pattern_repeat_freq
            anal_error_compare = (squeeze(mean(Err_norm_smoother(:,:,(i-1)*nsteps+1:i*nsteps,1), 1)) ...
                - squeeze(mean(Err_norm_anal_1(:,i,:,1:end-1,1), 1))) ...
                ./ squeeze(mean(Err_norm_anal_1(:,i,:,1:end-1,1), 1));
        else
            anal_error_compare = (squeeze(mean(Err_norm_smoother(:,:,(i-1)*nsteps+1:i*nsteps+1,1), 1)) ...
                - squeeze(mean(Err_norm_anal_1(:,i,:,:,1), 1))) ...
                ./ squeeze(mean(Err_norm_anal_1(:,i,:,:,1), 1));
        end
        plot((basetime:h:basetime+nsteps*h), bg_error_compare, 'Color', '#0072BD')
        if i ~= n_cycles_per_smoother*ob_pattern_repeat_freq
            plot((basetime:h:basetime+(nsteps-1)*h), anal_error_compare, 'r')
        else
            plot((basetime:h:basetime+nsteps*h), anal_error_compare, 'r')
        end
        xline(basetime + nsteps * h,'HandleVisibility','Off');
        yline(0,'HandleVisibility','Off')
    end
    title(sprintf('Relative change in average error norm - Total'))
    ylabel('Relative change')
    xlabel('Time within the window')
    hold off
    
    % Atmospheric norm
    subplot(3,1,2)
    hold on
    for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
        basetime = (i-1) * nsteps * h;
        bg_error_compare = (squeeze(mean(Err_norm_bg_2(:,:,i,:,2), 1)) ...
            - squeeze(mean(Err_norm_bg_1(:,i,:,:,2), 1))) ...
            ./ squeeze(mean(Err_norm_bg_1(:,i,:,:,2), 1));
        anal_error_compare = (squeeze(mean(Err_norm_anal_2(:,:,i,:,2), 1)) ...
            - squeeze(mean(Err_norm_anal_1(:,i,:,:,2), 1))) ...
            ./ squeeze(mean(Err_norm_anal_1(:,i,:,:,2), 1));
        plot((basetime:h:basetime+nsteps*h), bg_error_compare, 'Color', '#0072BD')
        plot((basetime:h:basetime+nsteps*h), anal_error_compare, 'r')
        xline(basetime + nsteps * h,'HandleVisibility','Off');
        yline(0,'HandleVisibility','Off')
    end
    title(sprintf('Relative change in average error norm - Atmosphere'))
    ylabel('Relative change')
    xlabel('Time within the window')
    hold off
    
    % Ocenaic norm
    subplot(3,1,3)
    hold on
    for i = 1:n_cycles_per_smoother*ob_pattern_repeat_freq
        basetime = (i-1) * nsteps * h;
        bg_error_compare = (squeeze(mean(Err_norm_bg_2(:,:,i,:,3), 1)) ...
            - squeeze(mean(Err_norm_bg_1(:,i,:,:,3), 1))) ...
            ./ squeeze(mean(Err_norm_bg_1(:,i,:,:,3), 1));
        anal_error_compare = (squeeze(mean(Err_norm_smoother(:,:,(i-1)*nsteps+1:i*nsteps+1,2), 1)) ...
            - squeeze(mean(Err_norm_anal_1(:,i,:,:,3), 1))) ...
            ./ squeeze(mean(Err_norm_anal_1(:,i,:,:,3), 1));
        plot((basetime:h:basetime+nsteps*h), bg_error_compare, 'Color', '#0072BD')
        plot((basetime:h:basetime+nsteps*h), anal_error_compare, 'r')
        xline(basetime + nsteps * h,'HandleVisibility','Off');
        yline(0,'HandleVisibility','Off')
    end
    title(sprintf('Relative change in average error norm - Ocean'))
    ylabel('Relative change')
    xlabel('Time within the window')
    hold off
end

if save_all_figures == 1
    addpath('C:\GitHub\PDRA-project\coupled L96\savepicspackage')
    mkdir(dirname);
    figHandles = findall(0,'Type','figure');
    for i = 1:numel(figHandles)
        fig_num = figHandles(i).Number;
        fn = strcat(dirname,num2str(fig_num),'-',num2str(i_trial),'-',...
            num2str(i_ob_pattern_repeats));  %in this example, we'll save to a temp directory.
        export_fig(fn,'-png',figHandles(i))
    end
    parameters_save.var_b = [var_atmos_bg,var_ocean_bg];
    parameters_save.var_o = var_ob;
    parameters_save.H = H;
    parameters_save.B = B;
    save(strcat(dirname,'parameters'),'parameters_save')
end


