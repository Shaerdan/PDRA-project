function [z_out] = f_model_l96c(z,no,na,alph,gamma,Fx,Fy)
%    "The Coupled Lorenz 1996 model."
%    z_out is the rhs vector of the coupled l96 model, computed at the
%    state z.
z_out = zeros(no+na,1);
y_out = zeros(no,1);
x_out = zeros(no,1);
x = z(1:na);
y = z(na+1:end);
for j = 1:na
    ip1 = j+1; if ip1>na, ip1=ip1-na; end
    im1 = j-1; if im1<1, im1=im1+na; end
    im2 = j-2; if im2<1, im2=im2+na; end
    x_out(j) = -x(im1)*(x(im2)-x(ip1))-x(j)...
        +Fx -  gamma *  y(j);
    y_out(j) = -alph*alph*y(im1)*(y(im2) - y(ip1)) -...
        alph*y(j) + Fy + gamma * x(j);
end
z_out = [x_out; y_out];

end

