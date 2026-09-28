linux=1;
digits(64);
miu0=pi*4e-7;
epsilon0=8.85e-12;
freq=5e13*(1:200)/200;
sigma=max(freq)/200;
weight_par=zeros(200,nq);
weight_perp=zeros(200,nq);
weight_a=zeros(200,nq);
weight_b=zeros(200,nq);
weight_c=zeros(200,nq);
fidin1=fopen( 'neig.txt','r');
fidout=fopen('Edispersion.txt','w');
fidout2=fopen('ovlp.txt','w');       % 与声子的内积
fidout3=fopen('dispeigvec.txt','w'); % Php的本征矢
fprintf(fidout3,'%d\n',nband);
fprintf(fidout3,'%d\n',nq);
% 看有几种原子
itype=1;
type_mass(itype)=mass_au(1);
for i=2:natom
    newtype=1;
    for j=1:i-1
        if abs(mass_au(i)-mass_au(j))<0.5
            newtype=0;
        end
    end
    if newtype==1
            itype=itype+1;
            type_mass(itype)=mass_au(i);
    end
end
ntype=itype
type_mass
for i=1:natom
    for itype=1:ntype
        if abs(mass_au(i)-type_mass(itype))<0.5
            type(i)=itype;
        end
    end
end
type
num=100000;
iq(1:num)=0;
diste(1:num)=0;
eigenw(1:num)=0;
fidin2=fopen('dispersion.txt','r');
i=1;
iq(i)=fscanf(fidin2,'%d',1);
diste(i)=fscanf(fidin2,'%g',1);
eigenw(i)=fscanf(fidin2,'%g',1);
while ~feof(fidin2)
    i=i+1;
    iq(i)=fscanf(fidin2,'%d',1);
    diste(i)=fscanf(fidin2,'%g',1);
    eigenw(i)=fscanf(fidin2,'%g',1);
    if i>num
        disp('num should be > the number of lines of the data file')
        stop
    end
end
fclose(fidin2);
num=i;
for i=1:num
    qnorm(iq(i))=norm(qp(iq(i),:));
    neig=fscanf(fidin1,'%d',1);
    degen=1;
       ow=eigenw(i)*2*pi;
        if ow>=0
            w=ow^2;
        else
            w=-ow^2;
        end
        qcross=[-qp(iq(i),2)^2-qp(iq(i),3)^2,qp(iq(i),1)*qp(iq(i),2),qp(iq(i),1)*qp(iq(i),3);qp(iq(i),1)*qp(iq(i),2),-qp(iq(i),1)^2-qp(iq(i),3)^2,qp(iq(i),2)*qp(iq(i),3);qp(iq(i),1)*qp(iq(i),3),qp(iq(i),2)*qp(iq(i),3),-qp(iq(i),1)^2-qp(iq(i),2)^2];
        DI0=reshape(DM(iq(i),:,:),nband,nband);
        DI=DI0-w*eye(nband);
        % DI是否可逆
        scale=1;
        for iband=1:nband
            scale=scale*norm(DI(iband,:))/1e28;
        end
        if 1%det(DI/1e28)>scale/1e200   % DI 可逆
            A=1.602e-19^2/Vol*Zm*(DI\(Zm.'));
            miuw2=miu0*w;
            emiuw2=epsilon0*miuw2;
            B=miuw2*A+emiuw2*eps0+qcross;
            Efield2=0;Efield3=0;
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
            [~,ii]=min(abs(Efield));
            if ii==1
                temp=[1 0 0];
            elseif ii==2
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
                if degen==1
                    Efield2=temp;
                    degen=degen+1;   % 与前两个方向垂直的电场若能使eq. (10)左边也很小，简并度再加1
                elseif degen==2
                    Efield3=temp;
                    degen=degen+1;
                end
            end
            if degen>2
                %disp('Warning: triple degenerate');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 至此求出了电场的两个本征态：Efield和Efield2，后者如果太小则不是。 %
            % 还可能有Efield3，具体看degen
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算原子振动
            dspl=1.602e-19*(DI\(Zm.'))*Efield.';   % 这是本征矢
            % 归一化
            norm1=norm(dspl);
            dspl=dspl/norm1;
            %==============将归一化的本征矢输出===============
            fprintf(fidout3,'%d  %g  %g  %g  %g\n',iq(i),qp(iq(i),1:3),ow/2/pi);
            for iband=1:nband
                fprintf(fidout3,'%g  %g  ',real(dspl(iband)),imag(dspl(iband))); 
            end
            fprintf(fidout3,'\n');% 每一行是一个本征矢
            fprintf(fidout3,'%g  %g  %g\n',Efield(1),Efield(2),Efield(3)); 
            %==================================================
            % 与band.yaml比较
            ovlp=ph_dspl*dspl;
            ovlp=abs(ovlp);
            %[~,i]=max(abs(ovlp)); 
            % 把dspl保存入dd然后调用plot_displacement.m
            iwp=1;nwp=1;
            dd(:,iwp)=dspl; % 定义iwp和dd只是为了迁就plot_displacement
            if iq(i)==7
                plot_displacement_BN % 画振动图
            else
            end
            if ow>1e11
                for j=1:nband
                    if ovlp(j)>0.1
                        fprintf(fidout2,'%g  %g  %g  %g\n',norm(qp(iq(i),:)),ow/2/pi,j,ovlp(j));
                    end
                end
            end
            if norm(Efield2)>1e-15
                dspl=1.602e-19*(DI\(Zm.'))*Efield2.';
                norm2=norm(dspl);
                dspl=dspl/norm2;
                %==============输出另一个可能的本征态===============
                fprintf(fidout3,'%d  %g  %g  %g  %g\n',iq(i),qp(iq(i),1:3),ow/2/pi);
                for iband=1:nband
                    fprintf(fidout3,'%g  %g  ',real(dspl(iband)),imag(dspl(iband))); 
                end
                fprintf(fidout3,'\n');% 每一行是一个本征矢
                fprintf(fidout3,'%g  %g  %g\n',Efield2(1),Efield2(2),Efield2(3)); 
                %==================================================
                ovlp=ph_dspl*dspl;
                ovlp=abs(ovlp);
                %[~,i]=max(abs(ovlp));
                % 准备把dspl保存入dd然后调用plot_displacement.m$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                iwp=2;nwp=2;
                dd(:,iwp)=dspl; % 定义iwp和dd只是为了迁就plot_displacement
                if iq(i)==7
                    plot_displacement_BN % 画振动图
                else
                end
                if ow>1e11
                    for j=1:nband
                        if ovlp(j)>0.1
                            fprintf(fidout2,'%g  %g  %g  %g\n',norm(qp(iq(i),:)),ow/2/pi,j,ovlp(j));
                        end
                    end
                end
            end
            if norm(Efield3)>1e-15
                dspl=1.602e-19*(DI\(Zm.'))*Efield3.';
                norm2=norm(dspl);
                dspl=dspl/norm2;
                %==============输出另一个可能的本征态===============
                fprintf(fidout3,'%d  %g  %g  %g  %g\n',iq(i),qp(iq(i),1:3),ow/2/pi);
                for iband=1:nband
                    fprintf(fidout3,'%g  %g  ',real(dspl(iband)),imag(dspl(iband))); 
                end
                fprintf(fidout3,'\n');% 每一行是一个本征矢
                fprintf(fidout3,'%g  %g  %g\n',Efield3(1),Efield3(2),Efield3(3)); 
                %==================================================
                ovlp=ph_dspl*dspl;
                ovlp=abs(ovlp);
                %[~,i]=max(abs(ovlp));
                % 准备把dspl保存入dd然后调用plot_displacement.m$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                iwp=2;nwp=2;
                dd(:,iwp)=dspl; % 定义iwp和dd只是为了迁就plot_displacement
                if iq(i)==7
                    plot_displacement_BN % 画振动图
                else
                end
                if ow>1e11
                    for j=1:nband
                        if ovlp(j)>0.1
                            fprintf(fidout2,'%g  %g  %g  %g\n',norm(qp(iq(i),:)),ow/2/pi,j,ovlp(j));
                        end
                    end
                end
            end
            Epar=dot(Efield,qp(iq(i),:))/norm(qp(iq(i),:));
            Eperp=cross(Efield,qp(iq(i),:))/norm(qp(iq(i),:));
            fprintf(fidout,'%g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n',norm(qp(iq(i),:)),ow/2/pi,abs(real(Efield(1:3)))/norm1,abs(imag(Efield(1:3)))/norm1,norm(Epar),norm(Eperp));
            if ow/2/pi>sigma*2
                for jjq=1:nq
                    for ifreq=1:200
                        iweight=gaussian_smearing(ow/2/pi,freq(ifreq),sigma)*gaussian_smearing(iq(i),jjq,0.1);%有些本征值没找全，所以有些亮点，此时需要把横向展展调大
                        weight_par(ifreq,jjq)=weight_par(ifreq,jjq)+iweight*norm(Epar)^2;
                        weight_perp(ifreq,jjq)=weight_perp(ifreq,jjq)+iweight*norm(Eperp)^2;
                        weight_a(ifreq,jjq)=weight_a(ifreq,jjq)+iweight*abs(Efield(1))^2;
                        weight_b(ifreq,jjq)=weight_b(ifreq,jjq)+iweight*abs(Efield(2))^2;
                        weight_c(ifreq,jjq)=weight_c(ifreq,jjq)+iweight*abs(Efield(3))^2;
                    end
                end
            end
            if norm(Efield2)>1e-15
                Epar=dot(Efield2,qp(iq(i),:))/norm(qp(iq(i),:));
                Eperp=cross(Efield2,qp(iq(i),:))/norm(qp(iq(i),:));
                fprintf(fidout,'%g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n',norm(qp(iq(i),:)),ow/2/pi,abs(real(Efield2(1:3)))/norm2,abs(imag(Efield2(1:3)))/norm2,norm(Epar),norm(Eperp));
                if ow/2/pi>sigma*2
                    for jjq=1:nq
                        for ifreq=1:200
                            iweight=gaussian_smearing(ow/2/pi,freq(ifreq),sigma)*gaussian_smearing(iq(i),jjq,0.1);%有些本征值没找全，所以有些亮点，此时需要把横向展展调大，或者算色散关系时要更大精度，找全本征值。
                            weight_par(ifreq,jjq)=weight_par(ifreq,jjq)+iweight*norm(Epar)^2;
                            weight_perp(ifreq,jjq)=weight_perp(ifreq,jjq)+iweight*norm(Eperp)^2;
                            weight_a(ifreq,jjq)=weight_a(ifreq,jjq)+iweight*abs(Efield2(1))^2;
                            weight_b(ifreq,jjq)=weight_b(ifreq,jjq)+iweight*abs(Efield2(2))^2;
                            weight_c(ifreq,jjq)=weight_c(ifreq,jjq)+iweight*abs(Efield2(3))^2;
                        end
                    end
                end
            end
        else % DI 不可逆
        end
    phononeig=real(sqrt(eig(DI0)))/2/pi;
    for iband=1:nband
        %subplot(2,2,3), plot(qnorm(iq(iiq)),phononeig(i),'ro');
        %hold on
    end
end
fprintf(fidout3,'-1');
fclose(fidout3);
figure
colormap('hot')
[X,Y] =  meshgrid(qnorm,freq);
subplot(2,3,1), pcolor(X/1e6,Y/1e12,sqrt(weight_par))
shading interp;
set(gca,'FontSize',18);
title('Parrallel')
xlabel('q (\mum^-^1)')
ylabel('Frequency (THz)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([0 temp(2)]);
subplot(2,3,2), pcolor(X/1e6,Y/1e12,sqrt(weight_perp))
shading interp;
set(gca,'FontSize',18);
title('Perpendicular')
xlabel('q (\mum^-^1)')
ylabel('Frequency (THz)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([0 temp(2)]);
subplot(2,3,4), pcolor(X/1e6,Y/1e12/0.03,sqrt(weight_a))
shading interp;
set(gca,'FontSize',18);
title('Direction a [001]')
xlabel('q (\mum^-^1)')
ylabel('Wave number (cm^-^1)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([450 950]);
subplot(2,3,5), pcolor(X/1e6,Y/1e12/0.03,sqrt(weight_b))
shading interp;
set(gca,'FontSize',18);
title('Direction b [100]')
xlabel('q (\mum^-^1)')
ylabel('Wave number (cm^-^1)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([450 950]);
subplot(2,3,6), pcolor(X/1e6,Y/1e12/0.03,sqrt(weight_c))
shading interp;
set(gca,'FontSize',18);
title('Direction c [010]')
xlabel('q (\mum^-^1)')
ylabel('Wave number (cm^-^1)')
temp=xlim;
xlim([0 temp(2)]);
temp=ylim;
ylim([450 950]);
colormap(hot)
fclose(fidout);
fclose(fidout2);
fclose(fidin1);


