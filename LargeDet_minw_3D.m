linux=1;
disp('allq 2024.5.20');
fidin=fopen('parameters.input','r');
ncpu=fscanf(fidin,'%d',1)
delete(gcp('nocreate'))
parpool(ncpu);
nw=fscanf(fidin,'%d',1)
minow=fscanf(fidin,'%g',1)
minow=minow*2*pi;
maxow=fscanf(fidin,'%g',1)   % Hz
maxow=maxow*2*pi;
fclose(fidin);
digits(64);
miu0=pi*4e-7;
epsilon0=8.85e-12;
fw(nw)=0;
ww(nw)=0;
zerop=zeros(nq,300);
nzero1=0*(1:nq);
nzero2=0*(1:nq);
fidout=fopen('dispersion.txt','w');
fid2d=fopen('dispersion2d.txt','w');
nsection=28;
for iq=1:nq
    qcross=[-qp(iq,2)^2-qp(iq,3)^2,qp(iq,1)*qp(iq,2),qp(iq,1)*qp(iq,3);qp(iq,1)*qp(iq,2),-qp(iq,1)^2-qp(iq,3)^2,qp(iq,2)*qp(iq,3);qp(iq,1)*qp(iq,3),qp(iq,2)*qp(iq,3),-qp(iq,1)^2-qp(iq,2)^2];
    DI0=reshape(DM(iq,:,:),nband,nband);
    for i=1:nw % par
        ow=minow+i*(maxow-minow)/nw;
        w1(i)=ow;
        % calculate fw
        if ow>=0
            w=ow^2;
        else
            w=-ow^2;
        end
        ww(i)=w;
        DI=DI0-w*eye(nband);
        DI(1:nband,nband+1:nband+3)=-1.602e-19*(Zm.');
        DI(nband+1:nband+3,1:nband)=miu0*w*1.602e-19/Vol*Zm;
        DI(nband+1:nband+3,nband+1:nband+3)=epsilon0*miu0*w*eps0+qcross;
        % 矩阵元素数量级相差太大，须平衡
        %DI(1:nband,:)=DI(1:nband,:)/max(max(abs(DI(1:nband,1:nband))));
        %DI(nband+1:nband+3,:)=DI(nband+1:nband+3,:)/max(max(abs(DI(nband+1:nband+3,1:nband))));
        factor=max(max(abs(DI(1:nband,nband+1:nband+3))));
        %DI(:,nband+1:nband+3)=DI(:,nband+1:nband+3)/factor;
        %[mintol,maxtol]=findtol_svd_ps(3,DI,nband);
        DI=DI0-w*eye(nband);
        fw(i)=det(DI);
    end
    if iq==2
        stop
    end
    ieig=0;
    for iw=2:nw-1
        if fw(iw)<fw(iw-1) && fw(iw)<fw(iw+1)    %  这是有非零解的必要条件之一
            ieig=ieig+1;
            eigenw(ieig)=w1(iw);
            iiw(ieig)=iw;
            % 因为频率太稀，在w1(iw-1)至w1(iw+1)之间搜索最小的fw，存为error(ieig)
            Lw=w1(iw-1); Lf=fw(iw-1);
            Rw=w1(iw+1); Rf=fw(iw+1);
            for iter=1:4
                parfor k=2:nsection+1%par
                    tw(k)=Lw+(Rw-Lw)/(nsection+1)*k;
                    w2=tw(k)^2;
                    DI=DI0-w2*eye(nband);
                    DI(1:nband,nband+1:nband+3)=-1.602e-19*(Zm.');
                    DI(nband+1:nband+3,1:nband)=miu0*w2*1.602e-19/Vol*Zm;
                    DI(nband+1:nband+3,nband+1:nband+3)=epsilon0*miu0*w2*eps0+qcross;
                    % 矩阵元素数量级相差太大，须平衡
                    DI(1:nband,:)=DI(1:nband,:)/max(max(abs(DI(1:nband,1:nband))));
                    DI(nband+1:nband+3,:)=DI(nband+1:nband+3,:)/max(max(abs(DI(nband+1:nband+3,1:nband))));
                    factor=max(max(abs(DI(1:nband,nband+1:nband+3))));
                    DI(:,nband+1:nband+3)=DI(:,nband+1:nband+3)/factor;
                    [mintol,maxtol]=findtol_svd_ps(3,DI,nband);
                    tfw(k)=mintol;
                end
                tw(1)=Lw;      tfw(1)=Lf;
                tw(nsection+2)=Rw; tfw(nsection+2)=Rf;
                [~,i]=min(tfw);
                if i==1
	                Rw=tw(2);
	                Rf=tfw(2);
	            elseif i==nsection+2
	                Lw=tw(nsection+1);
	                Lf=tfw(nsection+1);
	            elseif tfw(i-1)<tfw(i+1)
	                Lw=tw(i-1);
	                Lf=tfw(i-1);
	                Rw=tw(i);
	                Rf=tfw(i);
	            else
	                Lw=tw(i);
	                Lf=tfw(i);
	                Rw=tw(i+1);
	                Rf=tfw(i+1);
                end
            end
            error(ieig)=(Lf+Rf)/2;
            eigenw_precise(ieig)=(Lw+Rw)/2;
        end
    end
    if ieig==0
        disp('Warning: May be DI is not balanced');stop
    end
    neig=ieig;
    diste=sqrt(qp(iq,1)^2+qp(iq,2)^2+qp(iq,3)^2);
    for ieig=1:neig
        fprintf(fid,'\n');
        fprintf(fidout,'%d  %g  %.15g\n',iq,diste,eigenw_precise(ieig)/2/pi);
        fprintf(fid2d,'%g  %g  %g  %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw_precise(ieig)/2/pi);
    end
end
fprintf(fid2d,'-1.0  -1.0  -1,0  -1.0'); % Indicate the end of file
fclose(fidout);
fclose(fid2d);
fclose(fidout);