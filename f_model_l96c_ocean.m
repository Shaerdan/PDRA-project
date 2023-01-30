function [y_out] = f_model_l96c_ocean(y,x_c,no,alph,gamma,Fy)
%    "The Ocean model in Coupled L96 model with fixed atomsphere forcing x_c"
%     y_out is the rhs vector of such model.
%     y is the ocean state, x_c is the fixed atmosphere forcing.
y_out = zeros(no,1);
for j = 1:no
    ip1 = j+1; if ip1>no, ip1=ip1-no; end
    ip2 = j+2; if ip2>no, ip2=ip2-no; end
    im1 = j-1; if im1<1, im1=im1+no; end
    im2 = j-2; if im2<1, im2=im2+no; end
    y_out(j) = -alph*alph*y(im1)*(y(im2) - y(ip1)) -...
        alph*y(j) + Fy + gamma * x_c(j);
end

end

