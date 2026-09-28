clear
qp(17355,3)=0;
freq(17355)=0;
fid2d=fopen('dispersion2d.txt','r');
lbound=0;
ubound=40e12;
iq=1;
for i=1:3
    tq(i)=fscanf(fid2d,'%g',1);
end
qp(iq,1:3)=tq(:);
freq(iq)=fscanf(fid2d,'%g',1);
iq=iq+1;
while tq(1)>=0 && tq(2)>=0 
    for i=1:3
        tq(i)=fscanf(fid2d,'%g',1);
    end
    tf=fscanf(fid2d,'%g',1);
    if tf>lbound && tf<ubound
 %       if norm(tq(1:3)-qp(iq-1,:))<1e-6*norm(qp(iq-1,:))
 %           disp('Warning: this qp has two frequency.');
 %           if tf>freq(iq-1)  % 同一个q点，只取最大那个频率
 %               qp(iq-1,:)=tq(1:3);
 %               freq(iq-1)=tf;
 %           end
 %       else
            qp(iq,:)=tq(1:3);
            freq(iq)=tf;
            iq=iq+1;
 %       end
    end
end
freq=freq/1e12;
nq=iq-1;
mx=max(qp(:,1));
my=max(qp(:,2));
mz=max(qp(:,3));
qmax=min(my,mx);
ngrid=100;
% interp
clear x;
clear y;
theta=0/90*pi/2;
i=1;
for iq=1:nq
    angle=atan(qp(iq,1)/qp(iq,2));
    if(abs(angle-theta)<0.05)
        x(i)=sqrt(qp(iq,1)^2+(qp(iq,2)^2));
        y(i)=freq(iq);
        i=i+1;
    end
end
plot(x,y,'.')
xlabel('q (m^-^1)','FontSize',14)
ylabel('Frequency (THz)','FontSize',14)
title('fai=0 deg','FontSize',14)
axis([-inf inf 22 31])
