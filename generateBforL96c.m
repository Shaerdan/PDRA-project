function [Bc] = generateBforL96c(Nx,L,variance)
%   Construct climatological Bc from exponential kernal function
    xmin = 0; xmax = Nx-1; dx = 1;
    x_grid = xmin:dx:xmax+1;             
    Bc = zeros(Nx, Nx);
    L_grid = Nx*dx;

    for j =1:Nx
        for i = 1:Nx
            rad = (abs(x_grid(i)-x_grid(j)));
            rad = min(rad,L_grid-rad);        % Periodic domain
            Bc(i,j) = variance * exp(-rad/L);
        end
    end

end

