function [za2_f] = smoother_step(za_plot,assim_steps,s5_B_scaling,...
    Bo,Roinv,H,s5_smoother_loops,z_ob,n_cycles_per_smoother,l_integration_coupled_s5,s5_iterations,i_smooth_iteration,...
    h,nsteps,na,no,Fx,Fy,alph,gamma,ob_ix,i_ob_pattern_repeats,ob_pattern_repeat_freq,i_part_of_ob_pattern,l_lin_s5,...
    max_iterations,tolerance)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% za2_f assigment needs to be looked at!!!
%% Smoother
za_plot_2 = reshape(permute(za_plot(:,:,1:nsteps),[1,3,2]),[na+no,assim_steps]);
za_plot_2 = [za_plot_2 za_plot(:,end,end)];

Bo_smoother = s5_B_scaling * Bo;
Boinv_smoother = inv(Bo_smoother);

z_b = za_plot_2(:,1);

for i_count_smoother = 1:s5_smoother_loops
    if i_count_smoother == 1
        z_lin = za_plot_2;
    else
        z_lin = za2_f;
    end
    innov_o = z_ob(:,1:nsteps*n_cycles_per_smoother) - H*z_lin(:,2:nsteps*n_cycles_per_smoother+1);
    dX0_o=zeros(no,1);
    [dXo_anal2, JXoInner2, dJXoInner2, ito2, rel_grad_o2] = minimize_mod_crit_NKN(dX0_o,'calcfg_ocean_l96c',max_iterations,tolerance,...
        z_b,innov_o,z_lin,H,Boinv_smoother,Roinv,nsteps*n_cycles_per_smoother,h,na,no,Fx,Fy,alph,gamma,ob_ix(:,2));
    %                         assim_scheme,h,nsteps*n_cycles_per_smoother,beta,rho,sigma, ...
    %                         Omega,k,w_star,coupling_freq, z_b, B_smoother_inv, z_lin, 0, innov_o, var_ob5, ob_ix, l_fgat_s5);
    if l_integration_coupled_s5 && i_smooth_iteration == s5_iterations && i_count_smoother == s5_smoother_loops
        X_anal2 = z_lin(:,1) + [zeros(na,1);dXo_anal2];
        X_temp = X_anal2;
        for icycles = 1:n_cycles_per_smoother
            za2_f_icycles = l96c_rk2(X_temp,h,nsteps,na,no,Fx,Fy,alph,gamma);
            za2_f(na+1:na+no,(icycles-1)*nsteps+1:icycles*nsteps) = za2_f_icycles(na+1:end,1:nsteps);
            if icycles == n_cycles_per_smoother
                za2_f(na+1:na+no,end) = za2_f_icycles(na+1:na+no,nsteps+1);
                za2_f(1:na,:) = za_plot_2(1:na,:);
            end
            X_temp(1:na) = za_plot_2(1:na,icycles*nsteps+1); %
            X_temp(na+1:na+no) = za2_f(na+1:na+no,nsteps+1);
        end
    else
        X_anal2 = z_lin(:,1) + [zeros(na,1);dXo_anal2];
        [ya2_f] = l96c_rk2_ocean(za_plot_2(1:na,1:end-1),X_anal2(na+1:end),...
            h,nsteps*n_cycles_per_smoother,no,Fx,Fy,alph,gamma);
        za2_f(1:na,:)=za_plot_2(1:na,:);
        za2_f(na+1:na+no,:)=ya2_f;
    end
end
% za2_f is the guess trajectory for the next iteration
% (i_smooth_iteration), and is the final analysis if this
% is the last iteration

% if l_lin_s5 == 0
%     Bg_err_smoother((i_ob_pattern_repeats-1) * ...
%     ob_pattern_repeat_freq + i_part_of_ob_pattern, :, :) = za_plot_2(na+1:end,:) - z(na+1:end,:);
% end

end


