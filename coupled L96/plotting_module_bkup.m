                    if i_cycles == n_cycles_per_smoother && assim_scheme == 5
                        disp('Fix the presmoothing plots for the smoother method!')
                        %                     figure(200 + i_ob_pattern_repeats)
                        %                     plot(za_chk(indx_show,:),'k-','DisplayName','Analysis Forecast'); hold on;
                        %                     plot(z(indx_show,:),'r-','DisplayName','Ground Truth'); hold on;
                        %                     plot(zb_f_chk(indx_show,:),'b-','DisplayName','Background Forecast'); hold on;
                        %                     plot(z_ob_chk(indx_show,:),'go','DisplayName','Observation');hold on;
                        %                     xlabel('Assimilation Steps')
                        %                     legend show
                        %                     figure(1200 + i_ob_pattern_repeats)
                        %                     semilogy(error_norm_bg,'b-','DisplayName','Background Trajectory Error Norm'); hold on;
                        %                     semilogy(error_norm_analysis,'k-','DisplayName','Analysis Trajectory Error Norm')
                        %                     legend show
                    elseif i_part_of_ob_pattern == ob_pattern_repeat_freq && assim_scheme == 4
                        z_ob_chk_store = [zeros(ntotal,1) z_ob_chk_store];
                        z_ob_chk_store(z_ob_chk_store == 0) = nan;
                        figure(200 + i_ob_pattern_repeats)
                        plot(za_chk_store(indx_show,:),'k-','DisplayName','Analysis Forecast'); hold on;
                        plot(z_store(indx_show,:),'r-','DisplayName','Ground Truth'); hold on;
                        plot(zb_f_chk_store(indx_show,:),'b-','DisplayName','Background Forecast'); hold on;
                        plot(z_ob_chk_store(indx_show,:),'go','DisplayName','Observation');hold on;
                        xlabel('Assimilation Steps')
                        legend show
                        figure(1200 + i_ob_pattern_repeats)
                        plot(error_norm_bg,'b-','DisplayName','Background Trajectory Error Norm'); hold on;
                        plot(error_norm_analysis,'k-','DisplayName','Analysis Trajectory Error Norm')
                        legend show
                    end