linux=1;
disp('allq 2024.5.20');
fidin=fopen('parameters.input','r');
ncpu=fscanf(fidin,'%d',1)
delete(gcp('nocreate'))
parpool(ncpu);
wstep=fscanf(fidin,'%d',1)
fclose(fidin);
digits(64);
miu0=pi*4e-7;
epsilon0=8.85e-12;
maxow=(5e13*2*pi)^1;
nw=wstep;
fw(nw)=0;
ww(nw)=0;
zerop=zeros(nq,300);
nzero1=0*(1:nq);
nzero2=0*(1:nq);
fidout=fopen('dispersion.txt','w');
fid2d=fopen('dispersion2d.txt','w');
for iq=1:nq
    qcross=[-qp(iq,2)^2-qp(iq,3)^2,qp(iq,1)*qp(iq,2),qp(iq,1)*qp(iq,3);qp(iq,1)*qp(iq,2),-qp(iq,1)^2-qp(iq,3)^2,qp(iq,2)*qp(iq,3);qp(iq,1)*qp(iq,3),qp(iq,2)*qp(iq,3),-qp(iq,1)^2-qp(iq,2)^2];
    DI0=reshape(DM(iq,:,:),nband,nband);
    parfor i=1:nw
        ow=i*maxow/nw;
        % calculate fw
        if ow>=0
            w=ow^2;
        else
            w=-ow^2;
        end
        ww(i)=w;
        DI=DI0-w*eye(nband);
        %DI=DI*sym('1');
        A=1.602e-19^2/Vol*Zm*(DI\(Zm.'));
        %A=double(A);
        miuw2=miu0*w;
        emiuw2=epsilon0*miuw2;
        B=miuw2*A+emiuw2*eps0+qcross;
        fw(i)=-B(1,3)*B(2,2)*B(3,1)+B(1,2)*B(2,3)*B(3,1)+B(1,3)*B(2,1)*B(3,2)-B(1,1)*B(2,3)*B(3,2)-B(1,2)*B(2,1)*B(3,3)+B(1,1)*B(2,2)*B(3,3);
    end
    fw=double(fw);
    if iq==2
%        stop
    end
    if iq==nq
        if linux==0
            subplot(1,2,1), plot(sqrt(ww)/2/pi,real(fw),'b*')
        end
        hold on
        if linux==0
            subplot(1,2,1), plot(sqrt(ww)/2/pi,imag(fw),'ro')
        end
        fidfw=fopen('fw.txt','w');
        for i=1:nw
            fprintf(fidfw,'%.10g  %g  %g\n',sqrt(ww(i))/2/pi,real(fw(i)),imag(fw(i)));
        end
        fclose(fidfw);
    end
    fw=real(fw);
    %    find zero points
    j=0;
    for i=1:nw-2
        if fw(i)*fw(i+1)<=0
            nzero1(iq)=nzero1(iq)+1;
            j=j+1;
            if abs(fw(i))<abs(fw(i+1))
                zerop(iq,j)=ww(i);
            else
                zerop(iq,j)=ww(i+1);
            end
            zerop(iq,j)=ww(i)-fw(i)*(ww(i+1)-ww(i))/(fw(i+1)-fw(i));
            if iq==4000
                disp(i);disp(ww(i:i+1));disp(fw(i:i+1));
            end
        else
            if abs(fw(i+2))>abs(fw(i+1))&&abs(fw(i))>abs(fw(i+1)) && fw(i)*fw(i+2)>0%这种情况要认真考虑
                bb=((fw(i)-fw(i+2))*(ww(i)^2-ww(i+1)^2)-(fw(i)-fw(i+1))*(ww(i)^2-ww(i+2)^2))/2/((fw(i)-fw(i+2))*(ww(i)-ww(i+1))-(fw(i)-fw(i+1))*(ww(i)-ww(i+2)));
                aa=(fw(i)-fw(i+1))/(ww(i)-2*bb+ww(i+1))/(ww(i)-ww(i+1));
                cc=fw(i)-aa*(ww(i)-bb)^2;
                if aa*cc<0
                    nzero2(iq)=nzero2(iq)+1;
                    if bb-sqrt(-cc/aa)>ww(i) && bb+sqrt(-cc/aa)<ww(i+2)
                        j=j+1;
                        zerop(iq,j)=bb;
			if iq==4
			    disp(i);disp(ww(i:i+2));disp(fw(i:i+2));
			    disp(aa);disp(bb);disp(cc);
			end
                    end
                else
                    if abs(cc)<max(abs(fw(i)),abs(fw(i+2)))/1000;
                        nzero2(iq)=nzero2(iq)+1;
                        j=j+1;
                        zerop(iq,j)=bb;
			if iq==4
                            disp(i);disp(ww(i:i+2));disp(fw(i:i+2));
                            disp(aa);disp(bb);disp(cc);
			end
                    end
                end
            end
        end
    end
    n0=j;
    % plot
    hold on
    if iq==1
        diste=0.0;
    else
        diste=diste+sqrt((qp(iq,1)-qp(iq-1,1))^2+(qp(iq,2)-qp(iq-1,2))^2+(qp(iq,3)-qp(iq-1,3))^2);
    end
    for j=1:n0
        if linux==0
	    if zerop(iq,j)>=0
                subplot(1,2,2), plot(diste,sqrt(zerop(iq,j))/2/pi,'ro');
            else
	        subplot(1,2,2), plot(diste,-sqrt(-zerop(iq,j))/2/pi,'ro');
	    end
        end
	if zerop(iq,j)>=0
            fprintf(fidout,'%d  %g  %.15g\n',iq,diste,sqrt(zerop(iq,j))/2/pi);
            fprintf(fid2d,'%g  %g  %g  %g\n',qp(iq,1),qp(iq,2),qp(iq,3),sqrt(zerop(iq,j))/2/pi);
	else
	    fprintf(fidout,'%d  %g  %.15g\n',iq,diste,-sqrt(-zerop(iq,j))/2/pi);
        fprintf(fid2d,'%g  %g  %g  %g\n',qp(iq,1),qp(iq,2),qp(iq,3),-sqrt(zerop(iq,j))/2/pi);
	end
    end
    drawnow;
end
nzero1
nzero2
nzero=nzero1+nzero2
hold off
fprintf(fid2d,'-1.0  -1.0  -1,0  -1.0'); % Indicate the end of file
fclose(fidout);
fclose(fid2d);
fidout=fopen('neig.txt','w');
for iq=1:nq
    fprintf(fidout,'%g  ',nzero(iq));
end
fclose(fidout);