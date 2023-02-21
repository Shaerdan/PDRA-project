

%%%%%  plotting of smoother trajectories and error norms
figure(800 + i_ob_pattern_repeats)
kk_sub_count = 1;
for i_sub_show = H_space_pattern(end/2+1:end)
    subplot(3,3,kk_sub_count)
    plot(za2_f(i_sub_show,:),'k-','DisplayName','PostSmoother Analysis Forecast'); hold on;
    plot(za_chk(i_sub_show,:),'b-','DisplayName','Presmoother Analysis Forecast'); hold on;
    plot(z(i_sub_show,:),'r-','DisplayName','True State');hold on;
    obs_chk = [0 z_ob(i_sub_show,:)];
    obs_chk(obs_chk == 0) = nan;
    plot(obs_chk,'go','DisplayName','Observation'); hold on;
    xlabel('Assimilation Steps')
    %                     legend show
    kk_sub_count = kk_sub_count + 1;
end
error_norm_postsmoother = vecnorm(za2_f(na+1:end,:) - z(na+1:end,:))/N_Obs_Num_Spatial_o;
error_norm_presmoother = vecnorm(za_chk(na+1:end,:) - z(na+1:end,:))/N_Obs_Num_Spatial_o;
figure(1400+ i_ob_pattern_repeats)
semilogy(error_norm_postsmoother,'k-','DisplayName','Postsmoother error norm'); hold on;
semilogy(error_norm_presmoother,'b-','DisplayName','Presmoother error norm');
xlabel('steps')
ylabel('Oceanic error norm')
legend show


% figure(400 + i_ob_pattern_repeats)
% kk_sub_count = 1;
% sub_row = 4;
% len_obs = length(H_space_pattern(end/2+1:end));
% for i_sub_show = H_space_pattern(1:end/2)
%     subplot(3,3,kk_sub_count)
%     plot(za2_f(i_sub_show,:),'k-','DisplayName','PostSmoother Analysis Forecast'); hold on;
%     plot(za_chk(i_sub_show,:),'b-','DisplayName','Presmoother Analysis Forecast'); hold on;
%     plot(z(i_sub_show,:),'r-','DisplayName','True State');hold on;
%     obs_chk = [0 z_ob(i_sub_show,:)];
%     obs_chk(obs_chk == 0) = nan;
%     plot(obs_chk,'go','DisplayName','Observation'); hold on;
%     xlabel('Assimilation Steps')
%     %                     legend show
%     kk_sub_count = kk_sub_count + 1;
% end
% error_norm_postsmoother = vecnorm(za2_f(1:na,:) - z(1:na,:))/N_Obs_Num_Spatial_a;
% error_norm_presmoother = vecnorm(za_chk(1:na,:) - z(1:na,:))/N_Obs_Num_Spatial_a;
% figure(1400+ i_ob_pattern_repeats)
% semilogy(error_norm_postsmoother,'k-','DisplayName','Postsmoother error norm'); hold on;
% semilogy(error_norm_presmoother,'b-','DisplayName','Presmoother error norm');
% xlabel('steps')
% ylabel('Atmospheric error norm')
% legend show


             za_chk(:,(i_cycles-1)*nsteps+1:(i_cycles-1)*nsteps+nsteps+1) = za_f;
                    z_ob_chk = [z_ob];
                    z_ob_chk(z_ob_chk == 0) = nan;
                    if assim_scheme == 4
                        indx_str = (i_part_of_ob_pattern-1)*nsteps+1;
                        indx_end = (i_part_of_ob_pattern-1)*nsteps + nsteps;
                        za_chk_store(:,indx_str:indx_end+1) = za_chk;
                        z_store(:,indx_str:indx_end+1) = z;
                        zb_f_chk_store(:,indx_str:indx_end+1) = zb_f_chk;
                        z_ob_chk_store(:,indx_str:indx_end) = z_ob_chk;
                        error_norm_analysis(indx_str:indx_end+1) = vecnorm(za_chk - z)/N_Obs_Num_Spatial;
                        error_norm_bg(indx_str:indx_end+1) = vecnorm(zb_f_chk - z)/N_Obs_Num_Spatial;
                    else
                        indx_str = (i_cycles-1)*nsteps+1;
                        indx_end = (i_cycles-1)*nsteps + nsteps;
                        error_norm_analysis(indx_str:indx_end) = vecnorm(za_chk(:,indx_str:indx_end)...
                            -z(:,indx_str:indx_end))/N_Obs_Num_Spatial;
                        error_norm_bg(indx_str:indx_end) = vecnorm(zb_f_chk(:,indx_str:indx_end)...
                            - z(:,indx_str:indx_end))/N_Obs_Num_Spatial;
                    end
                    % debug plotting za_f vs z;
                    indx_show = 10;
                    
                    
%             if assim_scheme == 4
%                 za_plot_pattern_repeats(:,:,i_part_of_ob_pattern,i_ob_pattern_repeats) = permute(za_plot,[1,3,2]);
%                 zb_plot_pattern_repeats(:,:,i_part_of_ob_pattern,i_ob_pattern_repeats) = permute(zb_plot,[1,3,2]);
%             end                    