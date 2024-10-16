warning('off');

function ph1=invert(zh1,wv2inv)
  ph1=zh1.*wv2inv;
endfunction

function [u,v]=caluv(ph,k,l,trunc)
  u=-real(ifft2(i*l.*ph.*trunc));
  v=real(ifft2(i*k.*ph.*trunc));
endfunction

function qdot=advect(q,u,v,k,l)
  qdot=i*k.*fft2(u.*q)+i*l.*fft2(v.*q);
endfunction

t=0;
tc=0;

psimax=[];
ts=[];
stat=[];

qh1=fft2(q1);

dqh1dt_p=0;
dqh2dt_p=0;
dt0=dt;dt1=0;

%frm=0;
while t<=tmax+dt/2

  q1=real(ifft2(qh1.*trunc));
  if(sum(isnan(q1(:)))>0)
    break;
  end
  ph1=invert(qh1,wv2inv);
  [u1,v1]=caluv(ph1,k,l,trunc);
%  f=real(ifft2(fc.*trunc));
%  fcorr=fft2(f.*sign(q1));

%  dqh1dt=-advect(q1,u1,v1,k,l)-r*qh1+fcorr;
  dqh1dt=-advect(q1,u1,v1,k,l)-r*qh1+fc;
  delth=dth0*randn(size(k));
  fc=fc.*exp(i*delth);
  
  if(rem(tc,tpl)==0)
    ts=[ts,t];
    stat=[stat,[mean(q1(:).^2);sqrt(mean(u1(:).^2+v1(:).^2))]];
    imagesc(x,y,q1);title(sprintf('q : t = %5.1f',t));
    axis('xy','square');
    frm=mkpng(frm);%drawnow;
%    plot(ts,stat(2,:));
%    drawnow();
  end

  qh1=frc.*filtr.*(qh1+dt0*dqh1dt+dt1*dqh1dt_p);

  dqh1dt_p=frc.*dqh1dt;

  if tc==0
    dt0=1.5*dt;dt1=-0.5*dt;
  end
  tc=tc+1;
  t=tc*dt;

end

%plot(ts,stat(2,:));
mean(stat(2,:))
%plot(stat(5,:),stat(6,:),'o')
%axis([-L/2,L/2,-L/2,L/2])
if(frm>0)
  mkpng()
end

%figure(2)
%qg1spect;
%figure(1)
