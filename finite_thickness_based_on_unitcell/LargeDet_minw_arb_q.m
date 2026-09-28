fidin=fopen('parameters.input','r');
thick=fscanf(fidin,'%g',1)
ncpu=fscanf(fidin,'%d',1)
parpool(ncpu);
qxi=fscanf(fidin,'%g',1)
qx=qp(iq,1)+1j*qxi;
qy=qp(iq,2);
c=3e8;
miu0=pi*4e-7;
epsilon0=8.85e-12;
minlamda=fscanf(fidin,'%g',1)
lamdastep=fscanf(fidin,'%g',1)
maxlamda=fscanf(fidin,'%g',1)
nlam=(maxlamda-minlamda)/lamdastep
minqzr=fscanf(fidin,'%g',1)
qzrstep=fscanf(fidin,'%g',1)
maxqzr=fscanf(fidin,'%g',1)
nqzr=(maxqzr-minqzr)/qzrstep+1
nw=fscanf(fidin,'%d',1)
disp('frequency range read in unit Hz');
minw=fscanf(fidin,'%g',1)
maxw=fscanf(fidin,'%g',1)
outputunit=fscanf(fidin,'%g',1)  % 1:THz, 2:cm-1, 3:meV
is_slab=fscanf(fidin,'%d',1)
niter=fscanf(fidin,'%d',1)
polarization=fscanf(fidin,'%d',1)
if polarization ==1 
    disp('p polarization considered');
elseif polarization ==2 
    disp('s polarization considered');
elseif polarization == 4
    disp('p+s polarization considered')
else
    disp('polarization parameter error');
    stop
end
if is_slab==0
    thick_ratio=1;
else
    thick_ratio=lattV(3,3)/thick;
end
nsection=28;
%minw=minw*0.03e12;
%maxw=maxw*0.03e12; % switch to THz
fclose(fidin);
wstep=(maxw-minw)/nw;
qzr=0;
qzi=0;
w(nw)=0;
mfrdiff(nw)=0;
eigenw(500)=0;
iiw(500)=0;
error(500)=0;
eigenw_precise(500)=0;
fid=fopen(['surface_dispersionq',num2str(iq,'%04d'),'.txt'],'w');
fid2d=fopen(['surface_dispersion2dq',num2str(iq,'%04d'),'.txt'],'w');
ndimM=nband;
for iqzr=1:nqzr  
        iqzr
        qzr=minqzr+qzrstep*(iqzr-1);
    for ilam=1:nlam
        if mod(ilam,100)==0
            disp(datetime);
            disp(ilam);
        end
        qzi=minlamda+lamdastep*ilam;
        qz=qzr-1j*qzi;
        %qx=sqrt(qp(iq,1)^2+qp(iq,2)^2)+1j*qxi;  %2023.1.8 debug
        qvec=[qx qy qz];
        DMq
        DMa=DM1;
        qvec=[qx qy -qz];
        DMq
        DMb=DM1;
        parfor iw=1:nw  %%%%%%%%%%%%% parfor
            w(iw)=(minw+wstep*(iw-1))*2*pi;
            w2=w(iw)^2;
            if polarization ==1
                [DI,colfactor]=makeDI_p(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio);
            elseif polarization == 2
                [DI,colfactor]=makeDI_s(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio);
            elseif polarization ==4
                if abs(qx/qy)<1e-6 || abs(qy/qx)<1e-6
                    disp('Error: q can not parallel to prinicpal axis, when polarzation == 4');
                    stop
                end
                [DI,colfactor]=makeDI_arb_q(w2,miu0,epsilon0,qx,qy,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio);
            else
                disp('polarization parameter error');
                stop
            end
            [mintol,maxtol]=findtol_svd_arb_q(polarization,DI,ndimM);
            %fw(iw)=(mintol+maxtol)/2;
            fw(iw)=mintol;
            %
        end
%        if iq==1
%            diste=0.0;
%        else
            diste=sqrt(qp(iq,1)^2+qp(iq,2)^2+qp(iq,3)^2);
%        end
        ieig=0;
        for iw=2:nw-1
            if fw(iw)<fw(iw-1) && fw(iw)<fw(iw+1)    %  这是有非零解的必要条件之一
                ieig=ieig+1;
                eigenw(ieig)=w(iw);
                iiw(ieig)=iw;
                % 因为频率太稀，在w(iw-1)至w(iw+1)之间搜索最小的fw，存为error(ieig)
                Lw=w(iw-1); Lf=fw(iw-1);
                Rw=w(iw+1); Rf=fw(iw+1);
                for iter=1:niter
                    parfor k=2:nsection+1%par
                        tw(k)=Lw+(Rw-Lw)/(nsection+1)*k;
                        w2=tw(k)^2;
                        if polarization ==1
                            [DI,colfactor]=makeDI_p(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio);
                        elseif polarization == 2
                            [DI,colfactor]=makeDI_s(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio);
                        elseif polarization ==4
                            [DI,colfactor]=makeDI_arb_q(w2,miu0,epsilon0,qx,qy,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio);
                        else
                            disp('polarization parameter error');
                            stop
                        end
                        [mintol,maxtol]=findtol_svd_arb_q(polarization,DI,ndimM);
		                %tfw(k)=(mintol+maxtol)/2;
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
            disp('Warning: May be DI is not balanced');
        end
        %eigenw=real(eigenw);
        neig=ieig;
        if polarization == 1 || polarization ==2
            for ieig=1:neig
                fprintf(fid,'\n');
                if outputunit==1 % THz
                    fprintf(fid,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g',iq,diste,iqzr,ilam,iiw(ieig),eigenw(ieig)/2/pi,eigenw_precise(ieig)/2/pi,qzr,qzi,error(ieig)); % 下一步要用pick3Dmin来找出error(iqzr,ilam,iiw)的局部极小值
                    fprintf(fid2d,'%g  %g  %g   %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw(ieig)/2/pi);
                elseif outputunit==2 % cm-1
                    fprintf(fid,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g',iq,diste,iqzr,ilam,iiw(ieig),eigenw(ieig)/2/pi/0.03e12,eigenw_precise(ieig)/2/pi/0.03e12,qzr,qzi,error(ieig)); % 下一步要用pick3Dmin来找出error(iqzr,ilam,iiw)的局部极小值
                    fprintf(fid2d,'%g  %g  %g   %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw(ieig)/2/pi/0.03e12);
                elseif outputunit==3    % meV
                    fprintf(fid,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g',iq,diste,iqzr,ilam,iiw(ieig),eigenw(ieig)/2/pi/0.03e12*0.124,eigenw_precise(ieig)/2/pi/0.03e12*0.124,qzr,qzi,error(ieig)); % 下一步要用pick3Dmin来找出error(iqzr,ilam,iiw)的局部极小值
                    fprintf(fid2d,'%g  %g  %g   %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw(ieig)/2/pi/0.03e12);
                else
                    disp('Error');
                    outputunit
                end
            end
        elseif polarization == 4
            for ieig=1:neig
                fprintf(fid,'\n');
                if outputunit==1 % THz
                    fprintf(fid,'%d  %g  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g',iq,qx,qy,iqzr,ilam,iiw(ieig),eigenw(ieig)/2/pi,eigenw_precise(ieig)/2/pi,qzr,qzi,error(ieig)); % 下一步要用pick3Dmin来找出error(iqzr,ilam,iiw)的局部极小值
                    fprintf(fid2d,'%g  %g  %g   %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw(ieig)/2/pi);
                elseif outputunit==2 % cm-1
                    fprintf(fid,'%d  %g  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g',iq,qx,qy,iqzr,ilam,iiw(ieig),eigenw(ieig)/2/pi/0.03e12,eigenw_precise(ieig)/2/pi/0.03e12,qzr,qzi,error(ieig)); % 下一步要用pick3Dmin来找出error(iqzr,ilam,iiw)的局部极小值
                    fprintf(fid2d,'%g  %g  %g   %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw(ieig)/2/pi/0.03e12);
                elseif outputunit==3    % meV
                    fprintf(fid,'%d  %g  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g',iq,qx,qy,iqzr,ilam,iiw(ieig),eigenw(ieig)/2/pi/0.03e12*0.124,eigenw_precise(ieig)/2/pi/0.03e12*0.124,qzr,qzi,error(ieig)); % 下一步要用pick3Dmin来找出error(iqzr,ilam,iiw)的局部极小值
                    fprintf(fid2d,'%g  %g  %g   %g\n',qp(iq,1),qp(iq,2),qp(iq,3),eigenw(ieig)/2/pi/0.03e12);
                else
                    disp('Error');
                    outputunit
                end
            end
        end
    end
end
fclose(fid);
linux=1;
if linux==1
    disp('Exit successufully!');
    iq
else
    plot(1,1);
end
%surf(lamda,w/2/pi,log(tol))
%view(0,90)
%shading interp
