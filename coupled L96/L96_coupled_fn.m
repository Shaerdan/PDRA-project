function [z_out] = L96_coupled_fn(z0,N,Fx,Fy,alpha,gamma)
%
% Integrate set of coupled L96 equations
%
x=zeros(N,2);
y=zeros(N,2);
x_out=zeros(N,1);
y_out=zeros(N,1);
%
x(:,1) = z0(1:N);
y(:,1) = z0(N+1:2*N);
%
for i=1:N
    ip1 = i+1; if ip1>N, ip1=ip1-N; end
    ip2 = i+2; if ip2>N, ip2=ip2-N; end
    im1 = i-1; if im1<1, im1=im1+N; end
    im2 = i-2; if im2<1, im2=im2+N; end
    x_out(i) = (-x(im2,1)*x(im1,1) + x(im1,1)*x(ip1,1) - x(i,1) + Fx  -  gamma *  y(i,1));
    y_out(i) = (-alpha*alpha*y(im2,1)*y(im1,1) + alpha*alpha*y(im1,1)*y(ip1,1) - alpha*y(i,1) + Fy + gamma * x(i,1));   
end
z_out = [x_out; y_out];