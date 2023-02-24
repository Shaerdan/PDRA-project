clc;
clear;
for ii = 1:1
    close all;
    n_repeat_no_cycle = 1;
    n_ob_pattern_repeats = 100;
    compare_or_standalone = [2,2]; % [1,1] for standalone weakly coupled 4d-var, 
                                   % [2,2] for standalone smoother, [1,2] for both and compare
    var_size_scaling = 1;                               
    var_bg = var_size_scaling*[0.1^2, 0.3^2];
    var_ratio_bgob = var_size_scaling*0.1;
    var_ob = var_ratio_bgob*[0.1^2, 0.3^2];
    space_skip_a = 1;
    space_skip_o = 2;
    ks = num2str([space_skip_a,space_skip_o]);
    ks = ks(ks ~= ' ');
    for update_method = 1:4
        tolerance = 1.0d-6;   % For solver
        max_iterations = 50; % For solver
        % Grad_Test = 1: turn on gradient tests for calcfg routines; = 0 turn off.
        save_all_figures = 1;
        if save_all_figures == 1
            dirname = strcat('C:\results\DA\21Feb2023\a',num2str(space_skip_a),...
            'o',num2str(space_skip_o),'-sigmaBsigmaR-',num2str(sqrt(var_bg)),'-bgobratio-',num2str(var_ratio_bgob),'\');
            dirname = dirname(dirname~=' ');
            mkdir(dirname);
            filename1 = strcat(dirname,'plot1a-',num2str(ii),'update_method',num2str(update_method),'.fig');
            filename2 = strcat(dirname,'plot1b-',num2str(ii),'update_method',num2str(update_method),'.fig');
        end
        nsteps = 4;
        h=0.0125d0;
        Fx=15;
        Fy=8;
        alph=0.5;
        gamma= 0.6;
        N = 40;
        na = N; no = N; ntotal = na + no;
        outer_loops = 1;     % number of outerloop for weakly coupled standard 4dvar
        s5_smoother_loops = 1;  % Number of outer loops for smoother step only
        
        H_space_pattern_atmos = 1:space_skip_a:na;
        H_space_pattern_ocean = 1:space_skip_o:na;
        
        H_diag_atmos = zeros(1,na); H_diag_ocean = zeros(1,na);
        H_diag_atmos(H_space_pattern_atmos) = 1; H_diag_ocean(H_space_pattern_ocean) = 1;
        H_atmos = diag(H_diag_atmos); H_ocean = diag(H_diag_ocean);
        % collapsing H by deleting zero only rows;
        H_atmos( all(~H_atmos,2), : ) = []; H_ocean( all(~H_ocean,2), : ) = [];
        H = blkdiag(H_atmos,H_ocean);
        N_Obs_Num_Spatial = size(H,1);
        N_Obs_Num_Spatial_a = size(H_atmos,1);
        N_Obs_Num_Spatial_o = size(H_ocean,1);
        R_atmos = var_ob(1)*eye(N_Obs_Num_Spatial_a,N_Obs_Num_Spatial_a); R_ocean = var_ob(2)*eye(N_Obs_Num_Spatial_o,N_Obs_Num_Spatial_o);
        R = blkdiag(R_atmos,R_ocean);
        Rinv = inv(R);
        Rainv = inv(R_atmos);
        Roinv = inv(R_ocean);
        
        number_of_samples = ntotal; % full sample size for (likely) nonsingular B
        l_SpCov_SOAR = 1; % 0 for sampled covariance B, 1 for SOAR
        L_atmos = 4; L_ocean = 8; % make these input variable
        
        [Bainv,Boinv,Ba,Bo,B,SD] = GetCovMatriceB(number_of_samples,h,10,na,no,Fx,Fy,alph,gamma,...
            l_SpCov_SOAR,L_atmos, L_ocean,var_bg(1), var_bg(1));
        disp(strcat('condition of B matrices = ',num2str([cond(Bo) cond(Ba)])))
        
        %% rk2 solver parameters and model parameters for the l96 coupled model
        for i_trial = compare_or_standalone(1):compare_or_standalone(2)
            rng(ii)  % Fix random seed
            % loop controls:
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
            l_plot_convergence = 0;
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

            %% Setting up parameters for the assimilation
            
            xvals=1:na; % atmosphere grid indeces
            yvals=1:no; % ocean grid indeces
            
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
            x0_init=sin(xvals/(na-1)*2*pi);
            y0_init=cos(5*yvals/(no-1)*2*pi);
            for i_repeat_one_cycle = 1: n_repeat_no_cycle
                if n_repeat_no_cycle ~= 1
                    i_model_tn = randperm(1000,10) + 200;
                    [z_chk] = l96c_rk2([x0_init';y0_init'],h,i_model_tn(end)*assim_steps,na,no,Fx,Fy,alph,gamma);
                    x0_init = z_chk(1:na,end)';
                    y0_init = z_chk(na+1:end,end)';
                end
                x_ob = (nsteps:nsteps:ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
                y_ob = (ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps:1:...
                    ob_pattern_repeat_freq*n_cycles_per_smoother*nsteps);
                
                for i_ob_pattern_repeats = 1:n_ob_pattern_repeats
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
                        z_t=[x0_t,y0_t]'; % truth at t_{(k-1)*assim_steps+1}, k is the repeated obs pattern
                        X_t= z_t;
                        
                        % Generate background
                        if (i_ob_pattern_repeats == 1 && i_part_of_ob_pattern == 1)
                            noise = randn(ntotal,1);
                            z_b = z_t + sqrtm(B) * noise(1:ntotal);
                        elseif assim_scheme == 5
                            z_b = za2_f(:,assim_steps+1);
                        elseif assim_scheme == 4
                            z_b = za_plot(:,end,end);
                        end
                        %                     disp(strcat('ii =',num2str(ii),' Pattern Repeats = ',num2str(i_ob_pattern_repeats),' i_trial =',num2str(i_trial)))
                        % initialising za_plot
                        za_plot=zeros(ntotal,n_cycles_per_smoother,nsteps+1);
                        ob_ix=zeros(assim_steps,2);
                        obs_noise=randn(ntotal,assim_steps);
                        z_ob=zeros(N_Obs_Num_Spatial,assim_steps);
                        z_ob_a = zeros(N_Obs_Num_Spatial_a,assim_steps);
                        z_ob_o = zeros(N_Obs_Num_Spatial_o,assim_steps);
                        x_ob_local = x_ob((x_ob > (i_part_of_ob_pattern-1)*assim_steps) & (x_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
                        y_ob_local = y_ob((y_ob > (i_part_of_ob_pattern-1)*assim_steps) & (y_ob <= i_part_of_ob_pattern*assim_steps)) - (i_part_of_ob_pattern-1)*assim_steps;
                        for i=x_ob_local
                            ob_ix(i,1) = 1;
                            z_ob_a(:,i) = H_atmos*(z(1:na,i+1) + sqrt(var_ob(1))*obs_noise(1:na,i));
                        end
                        for i=y_ob_local
                            ob_ix(i,2) = 1;
                            z_ob_o(:,i) =  H_ocean*(z(na+1:ntotal,i+1) + sqrt(var_ob(2))*obs_noise(na+1:ntotal,i));
                        end
                        z_ob = [z_ob_a; z_ob_o];
                        
                        
                        %
                        %% Outer loops for smoother
                        
                        for i_smooth_iteration=1:s5_iterations
                            for i_cycles = 1:n_cycles_per_smoother
                                if i_smooth_iteration == 1
                                    %% Background forecast
                                    disp('********* START OF NEW CYCLE *********')
                                    X_b=z_b;
                                    X_t=z_t;
                                    disp(strcat('ii =',num2str(ii),' Pattern Repeats = ',num2str(i_ob_pattern_repeats),' i_trial =',num2str(i_trial)))
                                    
                                    [zb_f] = l96c_rk2(z_b,h,nsteps,na,no,Fx,Fy,alph,gamma);
                                    xb_f = zb_f(1:na,:);
                                    yb_f = zb_f(na+1:ntotal,:);
                                    zb_plot(:,i_cycles,:) = zb_f;
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
                                    innov = z_ob(:,(i_cycles-1)*nsteps+1:i_cycles*nsteps) - H*z_lin(:,2:nsteps+1); % Assume no ob at time zero.

                                    
                                    %% weakly coupled standard 4dvar, minimisation routine:
                                    % Atmosphere
                                    dX0_a=zeros(na,1);
                                    [dXa_anal, JXaInner, dJXaInner, ita, rel_grad_a] = minimize_mod_crit_NKN(dX0_a,'calcfg_atmos_l96c',max_iterations,tolerance, ...
                                        zb_plot(:,i_cycles,1),innov(1:N_Obs_Num_Spatial_a,:),z_lin,H_atmos,...
                                        Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix(:,1));
                                    
                                    % Ocean
                                    dX0_o=zeros(no,1);
                                    [dXo_anal, JXoInner, dJXoInner, ito, rel_grad_o] = minimize_mod_crit_NKN(dX0_o,'calcfg_ocean_l96c',max_iterations,tolerance, ...
                                        zb_plot(:,i_cycles,1),innov(N_Obs_Num_Spatial_a+1:end,:),z_lin,H_ocean,...
                                        Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,0,ob_ix(:,2));
                                    format short
                                    IT_Outer_ITa_Ito = [i_count, ita, ito];
                                    format longe
                                    rel_grad_a;
                                    rel_grad_o;
                                    Cost_RG_Atmos = [JXaInner, rel_grad_a];
                                    Cost_RG_Ocean = [JXoInner, rel_grad_o];
                                    %%%%%%%%%%%%
                                    % store the atmosphere increment for one of
                                    % the smoother udpate method:
                                    dXa_icycle(:,i_cycles) = dXa_anal;
                                    % Update fields
                                    delta_z0 = [dXa_anal; dXo_anal];
                                    z0 = z0 + delta_z0;
                                    
                                    %%%%%%%%%%%
                                    format long
                                    X_anal = z0;                                    
                                end % end outerloop
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                disp(strcat('ii =',num2str(ii),' Pattern Repeats = ',num2str(i_ob_pattern_repeats),' i_trial =',num2str(i_trial)))
                                
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
                                    Bo,Roinv,H,H_ocean,s5_smoother_loops,z_ob,n_cycles_per_smoother,...
                                    h,nsteps,na,no,N_Obs_Num_Spatial_a,Fx,Fy,alph,gamma,ob_ix,l_lin_s5,...
                                    max_iterations,tolerance,update_method,dXa_icycle);
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
                                bg_norm = vecnorm(rdivide(squeeze(zb_plot(:,i,:)) - z(:,(i-1)*nsteps+1:i*nsteps+1), ntotal));
                                anal_norm = vecnorm(rdivide(squeeze(za_plot(:,i,:)) - z(:,(i-1)*nsteps+1:i*nsteps+1), ntotal));
                                if i_trial == 1
                                    Err_norm_bg_1(i_store_norm, i_part_of_ob_pattern, i, :, 1) = bg_norm;
                                    Err_norm_anal_1(i_store_norm, i_part_of_ob_pattern, i, :, 1) = anal_norm;
                                elseif i_trial == 2
                                    Err_norm_bg_2(i_store_norm, i_part_of_ob_pattern, i, :, 1) = bg_norm;
                                    Err_norm_anal_2(i_store_norm, i_part_of_ob_pattern, i, :, 1) = anal_norm;
                                end
                                
                                % Atmospheric norm
                                bg_norm = vecnorm(rdivide(squeeze(zb_plot(1:na,i,:)) - z(1:na,(i-1)*nsteps+1:i*nsteps+1), na));
                                anal_norm = vecnorm(rdivide(squeeze(za_plot(1:na,i,:)) - z(1:na,(i-1)*nsteps+1:i*nsteps+1), na));
                                if i_trial == 1
                                    Err_norm_bg_1(i_store_norm, i_part_of_ob_pattern, i, :, 2) = bg_norm;
                                    Err_norm_anal_1(i_store_norm, i_part_of_ob_pattern, i, :, 2) = anal_norm;
                                elseif i_trial == 2
                                    Err_norm_bg_2(i_store_norm, i_part_of_ob_pattern, i, :, 2) = bg_norm;
                                    Err_norm_anal_2(i_store_norm, i_part_of_ob_pattern, i, :, 2) = anal_norm;
                                end
                                
                                % Oceanic norm
                                bg_norm = vecnorm(rdivide(squeeze(zb_plot(na+1:ntotal,i,:)) - z(na+1:ntotal,(i-1)*nsteps+1:i*nsteps+1), no));
                                anal_norm = vecnorm(rdivide(squeeze(za_plot(na+1:ntotal,i,:)) - z(na+1:ntotal,(i-1)*nsteps+1:i*nsteps+1), no));
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
                                smoother_norm = vecnorm(rdivide(za2_f - z, ntotal));
                                Err_norm_smoother(i_store_norm, i_part_of_ob_pattern, :, 1) = smoother_norm;
                                
                                % Oceanic norm
                                smoother_norm = vecnorm(rdivide(za2_f(na+1:ntotal,:) - z(na+1:ntotal,:), no));
                                Err_norm_smoother(i_store_norm, i_part_of_ob_pattern, :, 2) = smoother_norm;
                            end
                        end % plotting module storage end
                        
                    end  % i_part_of_ob_pattern
                    disp('long window end')
                    if any(isnan(X_anal))
                        disp(strcat('error: divergence at i =',num2str(num2str(i_ob_pattern_repeats))))
                        break;
                    end
                end     % i_ob_pattern_repeats
                if any(isnan(X_anal))
                    disp(strcat('error: divergence at i =',num2str(num2str(i_ob_pattern_repeats))))
                    break;
                end
            end % one cycle repeat trials (no long cycling)
            
            %%%%%%%%%%% Plotting module
            if l_plot_avg_error_norm && i_repeat_one_cycle == n_repeat_no_cycle
                %% Plot error norms averaged over (smoother) cycles
                if i_trial == 1
                    plot1a = figure(50);
                elseif i_trial == 2
                    plot1b = figure(55+update_method);
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
        if exist('plot1a','var') == 1 && save_all_figures == 1
            saveas(plot1a,filename1);
        end
        if exist('plot1b','var') == 1 && save_all_figures == 1
            saveas(plot1b,filename2);
        end
        if any(isnan(X_anal))
            disp(strcat('error: divergence at i =',num2str(num2str(i_ob_pattern_repeats))))
            disp(strcat('condition of B',num2str([cond(Ba) cond(Bo)])))
            break;
        end
    end
end
