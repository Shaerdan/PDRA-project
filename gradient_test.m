function [] = gradient_test(dXa_test,dXo_test,zb,innov,...
    z_lin,H,ob_ix,Bainv,Rainv,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[fa1,ga1] = calcfg_atmos_l96c(dXa_test,zb,innov,...
    z_lin,H,Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
hg_a = ga1/norm(ga1,2);
[fo1,go1] = calcfg_ocean_l96c(dXo_test,zb,innov,...
    z_lin,H,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
hg_o = go1/norm(go1,2);
for itest = 1:14
    a(itest) = 10^(-itest);
    [fa2,ga2] = calcfg_atmos_l96c(dXa_test+a(itest)*hg_a,zb,innov,...
        z_lin,H,Bainv,Rainv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
    Phi_a_atmos(itest) = (fa2 - fa1)/(a(itest)*hg_a'*ga1);
    [fo2,go2] = calcfg_ocean_l96c(dXo_test+a(itest)*hg_o,zb,innov,...
        z_lin,H,Boinv,Roinv,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix);
    Phi_a_ocean(itest) = (fo2 - fo1)/(a(itest)*hg_o'*go1);
end
figure(150)
loglog(a,abs(Phi_a_atmos-1),'k','DisplayName','Gradient Test for Atmospheric Costs and Gradient')
xlabel('Perturbeation Scaling')
ylabel('$|\frac{||J(dXa+\gamma*h) - J(dXa)||}{||\gamma*h*dJ||}-1|$','interpreter','latex')
legend show

figure(151)
loglog(a,abs(Phi_a_ocean-1),'k','DisplayName','Gradient Test for Oceanic Costs and Gradient')
xlabel('Perturbeation Scaling')
ylabel('$|\frac{||J(dXa+\gamma*h) - J(dXa)||}{||\gamma*h*dJ||}-1|$','interpreter','latex')
legend show
end

