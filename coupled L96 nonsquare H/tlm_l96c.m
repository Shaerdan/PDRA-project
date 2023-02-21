function [F] = tlm_l96c(u,na,no,alph,gamma)
% TLM for l96 coupled
% Output: F is the full jacobian of coupled l96;

nz = na + no;
F_atom = zeros(na,na); F_ocean = zeros(no,no);
u_atom = u(1:na); u_ocean = u(na+1:nz);
% F_atom_ocean = -gamma*eye(na,na);
% F_ocean_atom = gamma*eye(no,no);
for j = 1:na
    ip1 = j+1; if ip1>na, ip1=ip1-na; end
    ip2 = j+2; if ip2>na, ip2=ip2-na; end
    im1 = j-1; if im1<1, im1=im1+na; end
    im2 = j-2; if im2<1, im2=im2+na; end
    F_atom(j,im2) = -u_atom(im1) ;
    F_atom(j,im1) = -(u_atom(im2) - u_atom(ip1));
    F_atom(j,j) = -1;
    F_atom(j,ip1) = u_atom(im1);
end
for j = 1:no
    ip1 = j+1; if ip1>no, ip1=ip1-no; end
    ip2 = j+2; if ip2>no, ip2=ip2-no; end
    im1 = j-1; if im1<1, im1=im1+na; end
    im2 = j-2; if im2<1, im2=im2+na; end    
    F_ocean(j,im2) = -(alph^2)*u_ocean(im1) ;
    F_ocean(j,im1) = -(alph^2)*(u_ocean(im2) - u_ocean(ip1));
    F_ocean(j,j) = -alph;
    F_ocean(j,ip1) = (alph^2)*u_ocean(im1);
end
F = blkdiag(F_atom,F_ocean);
% F(1:na,na+1:nz) = F_atom_ocean;
% F(na+1:nz,1:na) = F_ocean_atom;
end