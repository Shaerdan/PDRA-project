function [] = plot_convergence(JXa,dJXa,JXo,dJXo,i_cycles,i_part_of_ob_pattern,n_cycles_per_smoother,...
                                            l_plot_convergence)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
end

