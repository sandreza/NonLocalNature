global sh ns nx

%rand("seed",3703);
t=0;

dth=0.4;
r=0.1;
la=0.2;

nx=256;dt=1/8;
uamp=1;

cph=0;
%cph= -0.4;;

%kappa=0;
%kappa=5e-4;
%kappa=1e-4;
%kappa=1e-3;
%uamp=0;kappa=2e-3;

dx=4*pi/nx;
x=[0.5:nx]*dx;
y=[0.5:nx]'*dx;

[xg,yg]=meshgrid(x,y);

rfac=exp(-r*dt);

mx=round(2.3*uamp*dt/dx);
sh=ones(2*mx+1,1)*[1:nx];
ns=mx+1;
for n=-mx:mx
  sh(n+ns,:)=shift(sh(n+ns,:),n);
end
sh=sh';

thx=2*pi*rand(3,1);
thy=2*pi*rand(3,1);

function nv=flow(v0,th,x,dy,dt)
  vfac=sqrt(2/(1+0.56^2+0.40^2));
%  nv=round(v0*dt/dy*(sin(x+th(1))+cos(2*x+1.5*th(2))));
  nv=round(v0*dt/dy*vfac*(cos(0.5*x+th(1))+0.56*cos(x+th(2))+0.40*cos(1.5*x+th(3))));
%  ny=round(v0*dt/dy*(0.25*sin(3*x+th))+0.75*sin(4*x-th/2));
endfunction

function cn=rph(c,nv)
  global sh ns nx grd
  cn=zeros(nx,nx);
  for i=1:nx
%    tmp=c(:,i)+gfac*grd(:,nv(i)+ns);
%    cn(:,i)=tmp(sh(:,nv(i)+ns));
    cn(:,i)=c(sh(:,nv(i)+ns),i);
  end
endfunction

function cn=diffuse(c,k,delc)
  n=size(c,1);np1=n+1;np2=n+2;
  ce=[c(:,n),c,c(:,1)];
  ce=[ce(n,:)-delc;ce;ce(1,:)+delc];
  cn=(1-4*k)*c+k*(ce(2:np1,1:n)+ce(2:np1,3:np2)+ce(1:n,2:np1)+ce(3:np2,2:np1));
endfunction
