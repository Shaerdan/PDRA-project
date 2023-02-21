function [za_chk_store,z_store,zb_f_chk_store,z_ob_chk_store] = trajectory_plotting(indx_show,assim_scheme,nsteps,n_cycles_per_smoother,...
    i_ob_pattern_repeats,i_part_of_ob_pattern,i_cycles,z_ob_chk,za_chk,z,zb_f_chk)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if i_cycles == n_cycles_per_smoother && assim_scheme == 5
    figure(200 + i_ob_pattern_repeats)
    plot(za_chk(indx_show,:),'k-*','DisplayName','Analysis Forecast'); hold on;
    plot(z(indx_show,:),'r-','DisplayName','Ground Truth'); hold on;
    plot(zb_f_chk(indx_show,:),'b:','DisplayName','Background Forecast'); hold on;
    plot(z_ob_chk(indx_show,:),'g*','DisplayName','Observation');
    xlabel('Assimilation Steps')
    legend show
else
    figure(200 + i_ob_pattern_repeats)
    indx_str = (i_part_of_ob_pattern-1)*nsteps+1;
    indx_end = (i_part_of_ob_pattern-1)*nsteps + nsteps;
    za_chk_store(:,indx_str:indx_end) = za_chk;
    z_store(:,indx_str:indx_end+1) = z;
    zb_f_chk_store(:,indx_str:indx_end) = zb_f_chk;
    z_ob_chk_store(:,indx_str:indx_end) = z_ob_chk;
    plot(za_chk_store(indx_show,:),'k-*','DisplayName','Analysis Forecast'); hold on;
    plot(z_store(indx_show,:),'r-','DisplayName','Ground Truth'); hold on;
    plot(zb_f_chk_store(indx_show,:),'b:','DisplayName','Background Forecast'); hold on;
    plot(z_ob_chk_store(indx_show,:),'g*','DisplayName','Observation');
    xlabel('Assimilation Steps')
    legend show
end
end

