fid=fopen('clean_qpoints.txt','r');
nq=fscanf(fid,'%d',1);
fclose(fid);
clear
fid2d=fopen('dispersion2d.txt','r');
lbound=30.00e12;
ubound=30.31e12;
iq=1;
for i=1:3
    tq(i)=fscanf(fid2d,'%g',1);
end
qp(iq,1:3)=tq(:);
freq(iq)=fscanf(fid2d,'%g',1);
iq=iq+1;
while tq(1)>0 && tq(2)>0 
    for i=1:3
        tq(i)=fscanf(fid2d,'%g',1);
    end
    tf=fscanf(fid2d,'%g',1);
    if tf>lbound && tf<ubound
        if norm(tq(1:3)-qp(iq-1,:))<1e-6*norm(qp(iq-1,:))
            disp('Warning: this qp has two frequency.');
        else
            qp(iq,:)=tq(1:3);
            freq(iq)=tf;
            iq=iq+1;
        end
    end
end
freq=freq/1e12;
nq=iq-1;
mx=max(qp(:,1));
my=max(qp(:,2));
mz=max(qp(:,3));
ngrid=100;
% interp
for i=1:ngrid
    for j=1:ngrid
        qx(i,j)=(i-0.5)/ngrid*mx;
        qy(i,j)=(j-0.5)/ngrid*my;
        mdist=1e200;
        for iq=1:nq
            dist=(qp(iq,1)-qx(i,j))^2+(qp(iq,2)-qy(i,j))^2;
            if dist<mdist
                mdist=dist;
                fq(i,j)=freq(iq);
            end
        end
    end
end
% expand
for i=1:ngrid
    for j=1:ngrid
        xqx(ngrid+i,ngrid+j)=qx(i,j);
        xqy(ngrid+i,ngrid+j)=qy(i,j);
        xfq(ngrid+i,ngrid+j)=fq(i,j);
        xqx(ngrid+1-i,ngrid+j)=-qx(i,j);
        xqy(ngrid+1-i,ngrid+j)=qy(i,j);
        xfq(ngrid+1-i,ngrid+j)=fq(i,j);
        xqx(ngrid+i,ngrid+1-j)=qx(i,j);
        xqy(ngrid+i,ngrid+1-j)=-qy(i,j);
        xfq(ngrid+i,ngrid+1-j)=fq(i,j);
        xqx(ngrid+1-i,ngrid+1-j)=-qx(i,j);
        xqy(ngrid+1-i,ngrid+1-j)=-qy(i,j);
        xfq(ngrid+1-i,ngrid+1-j)=fq(i,j);
    end
end
surf(xqx,xqy,xfq);
view(0,90)
shading interp
caxis([lbound/1e12 ubound/1e12])
axis([-mx mx -my my])
xlabel('q_x (m^-^1)','FontSize',14)
ylabel('q_y (m^-^1)','FontSize',14)
set(gca,'FontSize',14)
colormap(jet)
c=colorbar;
c.Label.String = 'Frequency (THz)';
