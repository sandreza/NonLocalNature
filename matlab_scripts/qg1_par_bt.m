L=4*pi;
nx=128;ny=128;
L=8*pi;
nx=256;ny=256;
W=L;

%dt=1/256/2/2/2;
%dt=1/2^10;
%tpl=1/dt/4;
%tmax=100;

kappa=5*1e-6;
%hypo=5*2e-4;
r=0.0;
%kappa=1;
F1=0;%6.25;
%famp=3750;
%famp=1000;
%famp=750;
%famp=300;

dx=L/nx;
dy=W/ny;

k0x=2*pi/L;
k0y=k0x;
[k,l]=meshgrid([0:nx/2,-nx/2+1:-1]*k0x,[0:ny/2,-ny/2+1:-1]*k0y);

wv2=k.*k+l.*l;
wv2inv=-1./(wv2+F1);
wv2inv(1,1)=0;
wv2h=1./wv2;
wv2h(1,1)=0;

switch(nx)
    case 128
      cphi = 0.69*pi;
    case 256 
      cphi = 0.715*pi;
    case 512
      cphi = 0.735*pi;
    otherwise
      cphi = 0.65*pi;
end

wvx=sqrt((k*dx).^2+(l*dy).^2);
filtr=exp(-18*(wvx-cphi).^7).*(wvx>cphi)+(wvx<=cphi);
filtr(isnan(filtr))=1;

kmax2=((nx/2-1)*k0x).^2;
trunc=(wv2<kmax2);

[x,y]=meshgrid([1/2:1:nx]/nx*L-L/2,[1/2:1:ny]/ny*W-W/2);

function [xm,ym]=findmin(q,x,y)
  [m,j]=min(min(q));
  xm=x(1,j);
  [m,i]=min(min(q'));
  ym=y(i,1);
endfunction

q1=cos(3*x-2*y).*sin(5*y)+cos(3*x+2*y).*cos(0.5*x);
q1=0*q1;

% forcing
K=sqrt(wv2);
th=rand(size(k))*2*pi;
%fc=30000*exp(-3*(K-5).^2).*exp(i*th);
fc=famp*filtr.*exp(i*th);fc(1,1)=0;
%fc=250*filtr.*exp(i*th);fc(1,1)=0;
dth0=sqrt(2*dt);

frc=exp(-kappa*dt*wv2.^2-hypo*dt*wv2h.^2);

%figure(3)
%plot(sqrt(wv2),-log(frc)/dt,'.');
%%plot(sqrt(wv2),kappa*dt*wv2.^2+hypo*dt*wv2inv.^2,'.')
%figure(1)
