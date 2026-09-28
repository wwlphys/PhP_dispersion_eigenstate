linux=1;
digits(64);
miu0=pi*4e-7;
epsilon0=8.85e-12;
freq=5e13*(1:200)/200;
sigma=max(freq)/200;
weight_par=zeros(200);
weight_perp=zeros(200);
fidin1=fopen( 'neig.txt','r');
fidin2=fopen( 'dispersion.txt','r');
fidout=fopen('Edispersion.txt','w');
fidout2=fopen('ovlp.txt','w');
for iq=1:nq
    qnorm(iq)=norm(qp(iq,:));
    neig=fscanf(fidin1,'%d',1);
    if iq==5000
        fidtest=fopen('test.txt','w');
        DI0=reshape(DM(iq,:,:),nband,nband);
        for ow = 0:1.7e13/1000:1.7e13
            w=(ow*2*pi)^2;
            DI=DI0-w*eye(nband);
            fprintf(fidtest,'%g  %g  %g\n',ow,norm(DI),norm(inv(DI)));
        end
        fclose(fidtest);
        phononeig=real(sqrt(eig(DI0)))/2/pi;
        for i=1:nband
            plot(qnorm(iq),phononeig(i),'ro');
            hold on
        end
        stop
    end
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
        if iq==190
            %disp(norm(inv(DI)));
        end
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
            normB(1)=norm(B(1,:));
            normB(2)=norm(B(2,:));
            normB(3)=norm(B(3,:));
            if sum(normB)-max(normB)<1e-2*max(normB)  %若只有一个方程是独立的
                [~,i]=max(normB);
                CC=B(i,:);
                [~,i]=max(abs(CC));
                if i==1
                    Efield(1)=-CC(2)/CC(1);
                    Efield(2)=1;
                    Efield(3)=0;
                    Efield2(1)=-CC(3)/CC(1);
                    Efield2(2)=0;
                    Efield2(3)=1;
                elseif i==2
                    Efield(1)=1;
                    Efield(2)=-CC(1)/CC(2);
                    Efield(3)=0;
                    Efield2(1)=0;
                    Efield2(2)=-CC(3)/CC(2);
                    Efield2(3)=1;
                else
                    Efield(1)=1;
                    Efield(2)=0;
                    Efield(3)=-CC(1)/CC(3);
                    Efield2(1)=0;
                    Efield2(2)=1;
                    Efield2(3)=-CC(2)/CC(3);
                end
            else
                Efield2=0;
                %看哪两个方程最独立
                if min(normB)<1e-2*max(normB)
                    if normB(1)<=normB(2) && normB(1)<=normB(3)
                        C(1,:)=B(2,:);
                        C(2,:)=B(3,:);
                    elseif normB(2)<=normB(1) && normB(2)<=normB(3)
                        C(1,:)=B(1,:);
                        C(2,:)=B(3,:);
                    else
                        C(1,:)=B(1,:);
                        C(2,:)=B(2,:);
                    end
                else
                    col12=abs(dot(B(1,:),B(2,:)))/normB(1)/normB(2);
                    col13=abs(dot(B(1,:),B(3,:)))/normB(1)/normB(3);
                    col23=abs(dot(B(2,:),B(3,:)))/normB(2)/normB(3);
                    if col12<=col13&&col12<=col23
                        C(1,:)=B(1,:);
                        C(2,:)=B(2,:);
                    elseif col13<=col12&&col13<=col23
                        C(1,:)=B(1,:);
                        C(2,:)=B(3,:);
                    else
                        C(1,:)=B(2,:);
                        C(2,:)=B(3,:);
                    end
                end
                %解这两个最独立的方程
                tempc(1)=C(1,1)*C(2,3)-C(2,1)*C(1,3);
                tempc(2)=C(1,2)*C(2,1)-C(1,1)*C(2,2);
                tempc(3)=C(1,3)*C(2,2)-C(1,2)*C(2,3);
                if abs(tempc(1))<=abs(tempc(2)) && abs(tempc(3))<=abs(tempc(2))
                    Efield(3)=1;
                    Efield(2)=Efield(3)*tempc(1)/tempc(2);
                    Efield(1)=Efield(3)*tempc(3)/tempc(2);
                elseif abs(tempc(2))<=abs(tempc(3)) && abs(tempc(1))<=abs(tempc(3))
                    Efield(1)=1;
                    Efield(2)=Efield(1)*tempc(1)/tempc(3);
                    Efield(3)=Efield(1)*tempc(2)/tempc(3);
                else
                    Efield(2)=1;
                    Efield(1)=Efield(2)*tempc(3)/tempc(1);
                    Efield(3)=Efield(2)*tempc(2)/tempc(1);
                end
            end
            temp(1)=norm(Efield);
            Efield=Efield/temp(1);
            temp(1)=norm(Efield2);
            if temp(1)>1e-15;
                Efield2=Efield2/temp(1);
            end
            % 计算原子振动
            dspl=1.602e-19*(DI\(Zm.'))*Efield.';
            % 归一化
            temp(1)=norm(dspl);
            dspl=dspl/temp(1);
            %Efield=Efield/temp(1);
            for i=1:natom
                for j=1:3
                    %dspl(i*3+j-3)=dspl(i*3+j-3)/Mass(i);
                end
            end
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
                            iweight=gaussian_smearing(ow/2/pi,freq(ifreq),sigma)*gaussian_smearing(iq,iiq,3);%有些本征值没找全，所以有些亮点，此时需要把横向展展调大，或者算色散关系时要更大精度，找全本征值。
                            weight_par(ifreq,iiq)=weight_par(ifreq,iiq)+iweight*norm(Epar)^2;
                            weight_perp(ifreq,iiq)=weight_perp(ifreq,iiq)+iweight*norm(Eperp)^2;
                        end
                    end
                end
            end
        else % DI 不可逆
        end
        if norm(qp(iq,:))>1.4e6 && ow>12.5e12
            %stop
        end
    end
    phononeig=real(sqrt(eig(DI0)))/2/pi;
    for i=1:nband
        %subplot(2,2,3), plot(qnorm(iq),phononeig(i),'ro');
        %hold on
    end
end
figure
colormap('jet')
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


