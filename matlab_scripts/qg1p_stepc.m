warning('off');
%kappa_c=kappa;%1e-4;
frc_c=exp(-kappa_c*dt*wv2);

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
dt0=dt;dt1=0;
cf1=cf(1,:);

%frm=0;
while t<=tmax+dt/2

  q1=real(ifft2(qh1.*trunc));
  if(sum(isnan(q1(:)))>0)
    break;
  end
  ph1=invert(qh1,wv2inv);
  [u1,v1]=caluv(ph1,k,l,trunc);
  c=real(ifft2(ch.*trunc));

  dqh1dt=-advect(q1,u1,v1,k,l)-r*qh1+fc;
  dchdt=-advect(c,u1,v1,k,l)+cfh;
  
  delth=dth0*randn(size(k));
  fc=fc.*exp(i*delth);
  
  if(rem(tc,tpl)==0)
    ts=[ts,t];
    al=mean(mean(c).*cf1)/mean(cf1.*cf1);
    stat=[stat,[mean(q1(:).^2);sqrt(mean(u1(:).^2+v1(:).^2));al]];
    if(t>tstart)
        imagesc(x,y,cbar);title(sprintf('c : t = %g',t));
    else
        imagesc(x,y,c);title(sprintf('c : t = %g',t));
	colorbar()
    end	
    axis('xy','equal');
    frm=mkpng(frm);%drawnow;
  end

  qh1=frc.*filtr.*(qh1+dt0*dqh1dt+dt1*dqh1dt_p);
  ch=frc_c.*filtr.*(ch+dt0*dchdt+dt1*dchdt_p);
  dqh1dt_p=frc.*dqh1dt;
  dchdt_p=frc_c.*dchdt;

  if(t>tstart)
    cbar=1/(nbar+1)*(cbar*nbar+c);
    nbar = nbar+1;
  elseif (t==tstart)
    cbar=c;
    nbar=1;
  end
  if tc==0
    dt0=1.5*dt;dt1=-0.5*dt;
  end
  tc=tc+1;
  t=tc*dt;

end
plot(ts,stat(3,:));
return

plot(ts,stat(2,:));
mean(stat(2,:))
%plot(stat(5,:),stat(6,:),'o')
%axis([-L/2,L/2,-L/2,L/2])
if(frm>0)
  mkpng()
end