linux=1;
digits(64);
miu0=pi*4e-7;
epsilon0=8.85e-12;
freq=5e13*(1:200)/200;
sigma=max(freq)/200;
weight_par=zeros(200,nq);
weight_perp=zeros(200,nq);
fidin1=fopen( 'neig.txt','r');
fidin2=fopen( 'dispersion.txt','r');
fidout=fopen('Edispersion.txt','w');
fidout2=fopen('ovlp.txt','w');       % 与声子的内积
fidout3=fopen('dispeigvec.txt','w'); % Php的本征矢
fprintf(fidout3,'%d\n',nband);
fprintf(fidout3,'%d\n',nq);
for iq=1:nq
    qnorm(iq)=norm(qp(iq,:));
    neig=fscanf(fidin1,'%d',1);
    degen=1;
%    fprintf(fidout3,'%d\n',neig);
    for ieig=1:neig
       fscanf(fidin2,'%g',1);
       ow=fscanf(fidin2,'%g',1)*2*pi;
        if ow>=0
            w=ow^2;
        else
            w=-ow^2;
        end
        qcross=[-qp(iq,2)^2-qp(iq,3)^2,qp(iq,1)*qp(iq,2),qp(iq,1)*qp(iq,3);qp(iq,1)*qp(iq,2),-qp(iq,1)^2-qp(iq,3)^2,qp(iq,2)*qp(iq,3);qp(iq,1)*qp(iq,3),qp(iq,2)*qp(iq,3),-qp(iq,1)^2-qp(iq,2)^2];
        DI0=reshape(DM(iq,:,:),nband,nband);
        DI=DI0-w*eye(nband);
        % DI是否可逆
        scale=1;
        for i=1:nband
            scale=scale*norm(DI(i,:))/1e28;
        end
        if 1%det(DI/1e28)>scale/1e200   % DI 可逆
            A=1.602e-19^2/Vol*Zm*(DI\(Zm.'));
            miuw2=miu0*w;
            emiuw2=epsilon0*miuw2;
            B=miuw2*A+emiuw2*eps0+qcross;
            Efield2=0;
            mtarget=1e100;
            maxt=0;
            for t=0:pi/100:pi    % 哪个方向的电场能使eq.(10)左边最接近零？
                for f=0:2*pi/100:2*pi
                    target=abs(B(1,1)*sin(t)*cos(f)+B(1,2)*sin(t)*sin(f)+B(1,3)*cos(t))^2+abs(B(2,1)*sin(t)*cos(f)+B(2,2)*sin(t)*sin(f)+B(2,3)*cos(t))^2+abs(B(3,1)*sin(t)*cos(f)+B(3,2)*sin(t)*sin(f)+B(3,3)*cos(t))^2;
                    if target<mtarget
                        mt=t;
                        mf=f;
                        mtarget=target;
                    end
                    if target>maxt
                        maxt=target;
                    end
                end
            end
            Efield(1)=sin(mt)*cos(mf);
            Efield(2)=sin(mt)*sin(mf);
            Efield(3)=cos(mt);    % 这个电场能使eq.(10)左边最接近零
            [~,i]=min(abs(Efield));
            if i==1
                temp=[1 0 0];
            elseif i==2
                temp=[0 1 0];
            else
                temp=[0 0 1];
            end
            temp=cross(temp,Efield);
            temp=temp/norm(temp);
            target=abs(B(1,:)*temp.')^2+abs(B(2,:)*temp.')^2+abs(B(3,:)*temp.')^2;
            if target<maxt/1e8 || target<mtarget*1.05
                Efield2=temp;   % 与Efield垂直的电场若能使eq. (10)左边也很小，说明有两个本征态
                degen=2;
            else
                degen=1;
            end
            temp=cross(temp,Efield);
            target=abs(B(1,:)*temp.')^2+abs(B(2,:)*temp.')^2+abs(B(3,:)*temp.')^2;
            if target<maxt/1e8 || target<mtarget*1.05
                Efield2=temp;
                degen=degen+1;   % 与前两个方向垂直的电场若能使eq. (10)左边也很小，简并度再加1
            end
            if degen>2
                disp('Warning: triple degenerate');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 至此求出了电场的两个本征态：Efield和Efield2，后者如果太小则不是。 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算原子振动
            dspl=1.602e-19*(DI\(Zm.'))*Efield.';   % 这是本征矢
            % 归一化
            temp(1)=norm(dspl);
            dspl=dspl/temp(1);
            %==============将归一化的本征矢输出===============
            fprintf(fidout3,'%d  %g  %g  %g  %g\n',iq,qp(iq,1:3),ow/2/pi);
            for i=1:nband
                fprintf(fidout3,'%g  %g  ',real(dspl(i)),imag(dspl(i))); 
            end
            fprintf(fidout3,'\n');% 每一行是一个本征矢
            fprintf(fidout3,'%g  %g  %g\n',Efield(1),Efield(2),Efield(3)); 
            %==================================================
            % 与band.yaml比较
            ovlp=ph_dspl*dspl;
            [~,i]=max(abs(ovlp));
            if ow>1e11
                for j=1:nband
                    if ovlp(j)>0.1
                        fprintf(fidout2,'%g  %g  %g  %g\n',norm(qp(iq,:)),ow/2/pi,j,ovlp(j));
                    end
                end
            end
            if norm(Efield2)>1e-15
                dspl=1.602e-19*(DI\(Zm.'))*Efield2.';
                temp(1)=norm(dspl);
                dspl=dspl/temp(1);
                %==============输出另一个可能的本征态===============
                fprintf(fidout3,'%d  %g  %g  %g  %g\n',iq,qp(iq,1:3),ow/2/pi);
                for i=1:nband
                    fprintf(fidout3,'%g  %g  ',real(dspl(i)),imag(dspl(i))); 
                end
                fprintf(fidout3,'\n');% 每一行是一个本征矢
                fprintf(fidout3,'%g  %g  %g\n',Efield2(1),Efield2(2),Efield2(3)); 
                %==================================================
                ovlp=ph_dspl*dspl;
                [~,i]=max(abs(ovlp));
                if ow>1e11
                    for j=1:nband
                        if ovlp(j)>0.1
                            fprintf(fidout2,'%g  %g  %g  %g\n',norm(qp(iq,:)),ow/2/pi,j,ovlp(j));
                        end
                    end
                end
            end
            Epar=dot(Efield,qp(iq,:))/norm(qp(iq,:));
            Eperp=cross(Efield,qp(iq,:))/norm(qp(iq,:));
            fprintf(fidout,'%g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n',norm(qp(iq,:)),ow/2/pi,real(Efield(1:3)),imag(Efield(1:3)),norm(Epar),norm(Eperp));
            if ow/2/pi>sigma*2
                for iiq=1:nq
                    for ifreq=1:200
                        iweight=gaussian_smearing(ow/2/pi,freq(ifreq),sigma)*gaussian_smearing(iq,iiq,0.1);%有些本征值没找全，所以有些亮点，此时需要把横向展展调大
                        weight_par(ifreq,iiq)=weight_par(ifreq,iiq)+iweight*norm(Epar)^2;
                        weight_perp(ifreq,iiq)=weight_perp(ifreq,iiq)+iweight*norm(Eperp)^2;
                    end
                end
            end
            if norm(Efield2)>1e-15
                Epar=dot(Efield2,qp(iq,:))/norm(qp(iq,:));
                Eperp=cross(Efield2,qp(iq,:))/norm(qp(iq,:));
                fprintf(fidout,'%g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n',norm(qp(iq,:)),ow/2/pi,real(Efield2(1:3)),imag(Efield2(1:3)),norm(Epar),norm(Eperp));
                if ow/2/pi>sigma*2
                    for iiq=1:nq
                        for ifreq=1:200
                            iweight=gaussian_smearing(ow/2/pi,freq(ifreq),sigma)*gaussian_smearing(iq,iiq,0.1);%有些本征值没找全，所以有些亮点，此时需要把横向展展调大，或者算色散关系时要更大精度，找全本征值。
                            weight_par(ifreq,iiq)=weight_par(ifreq,iiq)+iweight*norm(Epar)^2;
                            weight_perp(ifreq,iiq)=weight_perp(ifreq,iiq)+iweight*norm(Eperp)^2;
                        end
                    end
                end
            end
        else % DI 不可逆
        end
    end
    phononeig=real(sqrt(eig(DI0)))/2/pi;
    for i=1:nband
        %subplot(2,2,3), plot(qnorm(iq),phononeig(i),'ro');
        %hold on
    end
end
fprintf(fidout3,'-1');
fclose(fidout3);
figure
colormap('hot')
[X,Y] =  meshgrid(qnorm,freq);
subplot(1,2,1), pcolor(X/1e6,Y/1e12,sqrt(weight_par))
shading interp;
set(gca,'FontSize',18);
title('Parrallel')
xlabel('q (\mum^-^1)')
ylabel('Frequency (THz)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([0 temp(2)]);
subplot(1,2,2), pcolor(X/1e6,Y/1e12,sqrt(weight_perp))
shading interp;
set(gca,'FontSize',18);
title('Perpendicular')
xlabel('q (\mum^-^1)')
ylabel('Frequency (THz)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([0 temp(2)]);
colormap(hot)
fclose(fidout);
fclose(fidout2);
fclose(fidin1);
fclose(fidin2);


