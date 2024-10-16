rphinit;
kappa=5e-4
frm=-1;

cav=0;

tmax=5000;tpl=10;
t=0;cnt=1;
%kappa=5e-4;
%rphinit;

c=zeros(nx,nx);
s=0.1*ones(nx,1)*sin(k0*x);
frm=-1;

while t<tmax
  u=flow(uamp,thy,y,dx,dt);
  v=flow(uamp,thx-cph*t,x,dx,dt);
  c=rph(rph(c,v)',u)'+s;
  c=diffuse(c,kappa*dt/dx^2,0);
  t=t+dt;
  thx=thx+dth*(rand(3,1)-0.5);
  thy=thy+dth*(rand(3,1)-0.5);
  if rem(t,tpl)==0
    if t<=1000
      imagesc(x,y,c);
      axis('xy','equal');titleb(sprintf('%5.3f',t));
      frm=mkpng(frm);
    else
      cav = cav*(cnt-1)/cnt+c/cnt;
      cnt += 1;
      imagesc(x,y,cav);
      axis('xy','equal');titleb(sprintf('av %5.3f',t));
      frm=mkpng(frm);
    end
  end
end

cb=mean(cav);bb=s(1,:);
alpha=sum(bb.*cb)/sum(bb.*bb)
kcoeff=1/alpha/k0-kappa
plot(x,cb,x,bb)
