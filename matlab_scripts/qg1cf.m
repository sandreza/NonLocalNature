ch=0*wv2;
dchdt_p=0;
cf=cos(k0*x);
cfh=fft2(cf);

%tmax=500;
qg1p_stepc

cm=mean(cbar);
cf1=cf(1,:);
al=mean(cm.*cf1)/mean(cf1.*cf1);
xa=x(1,:);
plot(xa,cm,xa,al*cf1)
keff=1/al/k0^2-kappa

# 0.5  0.3097
# 1.0  0.2514
# 1.5  0.2184
# 2.0  0.1912
# 2.5  0.1665
# 3.0  0.1441
