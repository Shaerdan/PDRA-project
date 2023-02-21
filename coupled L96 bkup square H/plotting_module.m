    %% plots of error norms
    if (l_plot_error_norm || l_plot_avg_error_norm)
        %% Store norm values for later plotting
        for i=1:n_cycles_per_smoother
            % Total norm
            bg_norm = vecnorm(rdivide(squeeze(zb_plot(:,i,:)) - z(:,(i-1)*nsteps+1:i*nsteps+1), SD'));
            anal_norm = vecnorm(rdivide(squeeze(za_plot(:,i,:)) - z(:,(i-1)*nsteps+1:i*nsteps+1), SD'));
            Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 1) = bg_norm;
            Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 1) = anal_norm;
            
            % Atmospheric norm
            bg_norm = vecnorm(rdivide(squeeze(zb_plot(1:na,i,:)) - z(1:na,(i-1)*nsteps+1:i*nsteps+1), SD(1:na)'));
            anal_norm = vecnorm(rdivide(squeeze(za_plot(1:na,i,:)) - z(1:na,(i-1)*nsteps+1:i*nsteps+1), SD(1:na)'));
            Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 2) = bg_norm;
            Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 2) = anal_norm;
            
            % Oceanic norm
            bg_norm = vecnorm(rdivide(squeeze(zb_plot(na+1:end,i,:)) - z(na+1:end,(i-1)*nsteps+1:i*nsteps+1), SD(na+1:end)'));
            anal_norm = vecnorm(rdivide(squeeze(za_plot(na+1:end,i,:)) - z(na+1:end,(i-1)*nsteps+1:i*nsteps+1), SD(na+1:end)'));
            Err_norm_bg(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 3) = bg_norm;
            Err_norm_anal(i_ob_pattern_repeats, i_part_of_ob_pattern, i, :, 3) = anal_norm;
        end
        
        if (assim_scheme == 5)
            % Total norm
            smoother_norm = vecnorm(rdivide(za2_f - z, SD'));
            Err_norm_smoother(i_ob_pattern_repeats, i_part_of_ob_pattern, :, 1) = smoother_norm;
            
            % Oceanic norm
            smoother_norm = vecnorm(rdivide(za2_f(na+1:end,:) - z(na+1:end,:), SD(na+1:end)'));
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
    if assim_scheme == 5
        ua2_f_store(:,l_marker+1:r_marker+1) = za2_f(:,:);
    end