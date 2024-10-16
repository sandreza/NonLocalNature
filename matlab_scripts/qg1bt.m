% set up turbulent state at equilibrium
famp=300;
%famp=500;
%famp=750;
%%famp=1500;
%famp=10;
%%famp=0.1;
%famp=1;
%famp=10;
%famp=50;

hypo=1e-3;
%hypo=1e-4;
%hypo=1e-1;
%hypo=1e-2;


dt=1/32;
dt=1/64;
dt=1/256;

tpl=1.0/dt;
tmax=2000;

qg1_par_bt;
%dth0=dth0*10;
q1=randn(size(q1));
q1=0*q1;
frm=-1;
%famp
qg1p_step;
famp
%plot(ts,stat(2,:))

# im(real(ifft2(ph1.*trunc)))

ke=0.5*stat(2,:).^2;
plot(ts,ke);
mean(ke(1001:end))

save -binary qg1bt.mat

return
qg1p_step;plot(ts,stat(2,:))
