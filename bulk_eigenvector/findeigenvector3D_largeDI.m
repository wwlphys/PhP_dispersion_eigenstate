linux=1;
outputunit=1;
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
fid3=fopen('Spinz_a.txt','w');
fid5=fopen('Spinx_a.txt','w');
fid7=fopen('Spiny_a.txt','w');
fidout=fopen('Edispersion.txt','w');
fidout2=fopen('ovlp.txt','w');       % 与声子的内积
fidout3=fopen('dispeigvec.txt','w'); % Php的本征矢
fidout4=fopen('phase.txt','w');      % 第一个原子x方向振动和三个电场的相位
fidout5=fopen('PPP_CrystalEnergy.txt','w');
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
% 定义一系列natom个右旋矢量和左旋矢量，每个矢量2*natom行
for i=1:natom
    for j=1:natom*2
        right(i,j)=0;
        left(i,j)=0;
        mass_right(i,j)=0;
        mass_left(i,j)=0;
    end
    right(i,2*i-1)= 1/sqrt(2);
    right(i,2*i  )=1j/sqrt(2);
    left(i,2*i-1) =  1/sqrt(2);
    left(i,2*i  ) =-1j/sqrt(2);
    mass_right(i,2*i-1)= Mass(i);
    mass_right(i,2*i  )=1j*Mass(i);
    mass_left(i,2*i-1) = Mass(i);
    mass_left(i,2*i  ) =-1j*Mass(i);
end
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
    ow=eigenw(i)*2*pi;
    if ow>=0
        w=ow^2;
    else
        w=-ow^2;
    end
    qcross=[-qp(iq(i),2)^2-qp(iq(i),3)^2,qp(iq(i),1)*qp(iq(i),2),qp(iq(i),1)*qp(iq(i),3);qp(iq(i),1)*qp(iq(i),2),-qp(iq(i),1)^2-qp(iq(i),3)^2,qp(iq(i),2)*qp(iq(i),3);qp(iq(i),1)*qp(iq(i),3),qp(iq(i),2)*qp(iq(i),3),-qp(iq(i),1)^2-qp(iq(i),2)^2];
    DI0=reshape(DM(iq(i),:,:),nband,nband);
    DI=DI0-w*eye(nband);
    DI(1:nband,nband+1:nband+3)=-1.602e-19*(Zm.');
    DI(nband+1:nband+3,1:nband)=miu0*w*1.602e-19/Vol*Zm;
    DI(nband+1:nband+3,nband+1:nband+3)=epsilon0*miu0*w*eps0+qcross;
    % 矩阵元素数量级相差太大，须平衡
    colfactor(1:nband+3)=1;  % 记下列变换乘的系数，用于后面计算本征矢
    DI(1:nband,:)=DI(1:nband,:)/max(max(abs(DI(1:nband,1:nband))));
    DI(nband+1:nband+3,:)=DI(nband+1:nband+3,:)/max(max(abs(DI(nband+1:nband+3,1:nband))));
    factor=max(max(abs(DI(1:nband,nband+1:nband+3))));
    DI(:,nband+1:nband+3)=DI(:,nband+1:nband+3)/factor;
    colfactor(nband+1:nband+3)=colfactor(nband+1:nband+3)/factor;
    nrow=nband+3;
    ncol=nband+3;
    minnormx=0;
    for i1=1:ncol % 假设本征矢的第i1个分量为1，求本征矢，存在vector中
        if i1==1
            MA=DI(:,2:ncol);
        elseif i1==ncol
            MA=DI(:,1:ncol-1);
        else
            MA=[DI(:,1:i1-1) DI(:,i1+1:ncol)]; % MA是nrow行，ncol-1列
        end
        MB=-DI(:,i1);
        x=pinv(MA)*MB; %x是ncol-1行，1列
        normx(i1)=norm(x);
        vector(:,i1)=[x(1:i1-1);1;x(i1:ncol-1)]; % vector是ncol行，ncol列
    end
    colinear(1:ncol)=0;
    clear vec;
    ideg=0;
    for i1=1:ncol % 上述本征矢中norm比较小的是对的，若两个本征矢正交，则有简并
        if normx(i1)<inf && colinear(i1)==0
            ideg=ideg+1;
            vec(:,ideg)=vector(:,i1);
        end
        for j1=i1+1:ncol
            if abs(vector(:,i1).'*vector(:,j1))/norm(vector(:,i1))/norm(vector(:,j1))>0.1
                colinear(j1)=1;
            end
        end
    end
    ndeg=ideg;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % 至此求出了ndeg个简并本征态：vec(:,1:ndeg),重新代回原方程，不能使方程成立的解扔掉，剩下的存在vec2(:,1:ndeg2)中
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear vec2
    vec2(ncol,ndeg)=0;
    ideg2=0;
    clear error;
    error(ndeg)=0;
    for ideg=1:ndeg
        error(ideg)=norm(DI*vec(:,ideg));
    end
    minerror=min(error);       % 输入min(error)看看，min(error)太大说明这个频率不是PhP，是纯声子。
    for ideg=1:ndeg
        if error(ideg)<minerror*2
            ideg2=ideg2+1;
            vec2(:,ideg2)=vec(:,ideg);
        end
    end
    ndeg2=ideg2;
    % makeDI里调整了数量级，这里调回去。
    for ideg2=1:ndeg2
        vec2(1:ncol,ideg2)=vec2(:,ideg2).*colfactor(:);
    end
     % 计算原子振动
    iwp=0;
    for ideg=1:ndeg2
        dspl=vec2(1:nband,ideg);   % 这是本征矢的原子振动a部分
        % 归一化
        norma=norm(dspl);
        dspl=dspl/norma;
        %==============将归一化的本征矢保存===============
        iwp=iwp+1;
        wp(iwp)=eigenw(i)*2*pi;  % 标准单位的圆频率
        dd(:,iwp)=dspl;  
        if iq(i)==7 && ideg == 2 || iq(i)==26
            %plot_displacement_vector_hBN % 画振动图
            plot_displacement_BN_qx
        else
        end
        %==============原子振动归一化后的电场强度=========
        Efield=vec2(nband+1:nband+3,ideg);
        Efield(1:3)=Efield(1:3)/norma;
        %==============将归一化的本征矢输出===============
        fprintf(fidout3,'%d  %g  %g  %g  %g\n',iq(i),qp(iq(i),1:3),ow/2/pi);
        for iband=1:nband
            fprintf(fidout3,'%g  %g  ',real(dspl(iband)),imag(dspl(iband))); 
        end
        fprintf(fidout3,'\n');% 每一行是一个本征矢
        fprintf(fidout3,'%g  %g  %g\n',Efield(1),Efield(2),Efield(3)); 
        %==============输出第一个原子X方向振动和三个电场的相位=========
        fprintf(fidout4,'%d  %g  %g  %g  %g  ',iq(i),qp(iq(i),1:3),ow/2/pi);
        temp1=angle(dspl(1));
        if temp1>2
            temp1=temp1-2*pi;
        end
        temp3=angle(dspl(3));
        if temp3>2
            temp3=temp3-2*pi;
        end
        temp=angle(Efield(1:3));
        for itemp=1:3
            if temp(itemp)>2
                temp(itemp)=temp(itemp)-2*pi;
            end
        end
        fprintf(fidout4,'%g  %g  %g  %g  %g\n',temp1,temp3,temp(1:3));
        %==============输出复的极化强度PPP、晶体势能和动能随时间变化========
        if ow/2/pi/0.03e12>700 && ow/2/pi/0.03e12<750
            fprintf(fidout5,'%d  %g  %g  %g  %g Hz  %g cm-1\n',iq(i),qp(iq(i),1:3),ow/2/pi,ow/2/pi/0.03e12);
            PPP(1:3)=0;
            dsplScale=1e-12*sqrt(Mass(1));  % 乘上这个数，免得原子位移大得离谱，原子位移应小于0.01埃，所以dspl应小于1e-12*根号M
            for iatom=1:natom
                for idd=1:3     %
                    for jdd=1:3
                        PPP(jdd)=PPP(jdd)+dspl((iatom-1)*3+idd)*dsplScale*Zm(jdd,(iatom-1)*3+idd); % Zm中除的根号M跟dspl含有的根号M相抵消
                    end
                end
            end
            PPP=PPP/Vol*1.602e-19; % abc是晶向
            fprintf(fidout5,'Pa = %f + %f*i = %f * exp(i*%f) unit:C/m^2\n',real(PPP(1)),imag(PPP(1)),abs(PPP(1)),angle(PPP(1)));
            fprintf(fidout5,'Pb = %f + %f*i = %f * exp(i*%f) unit:C/m^2\n',real(PPP(2)),imag(PPP(2)),abs(PPP(2)),angle(PPP(2)));
            fprintf(fidout5,'Pc = %f + %f*i = %f * exp(i*%f) unit:C/m^2\n',real(PPP(3)),imag(PPP(3)),abs(PPP(3)),angle(PPP(3)));
            angleax=77/180*pi; % xyz是实验室坐标系，angleax是a跟x的夹角，x轴是光栅的短轴，Y轴是光栅的长轴
            Pxyz(1)= PPP(1)*cos(angleax)+PPP(2)*sin(angleax);
            Pxyz(2)=-PPP(1)*sin(angleax)+PPP(2)*cos(angleax);
            Pxyz(3)= PPP(3);
            fprintf(fidout5,'Rota angle: %f deg\n',angleax*180/pi);
            fprintf(fidout5,'Px = %f + %f*i = %f * exp(i*%f) unit:C/m^2\n',real(Pxyz(1)),imag(Pxyz(1)),abs(Pxyz(1)),angle(Pxyz(1)));
            fprintf(fidout5,'Py = %f + %f*i = %f * exp(i*%f) unit:C/m^2\n',real(Pxyz(2)),imag(Pxyz(2)),abs(Pxyz(2)),angle(Pxyz(2)));
            fprintf(fidout5,'Pz = %f + %f*i = %f * exp(i*%f) unit:C/m^2\n',real(Pxyz(3)),imag(Pxyz(3)),abs(Pxyz(3)),angle(Pxyz(3)));
            %输出随时间变化的晶体势能密度+动能密度
            fprintf(fidout5,'Time(s) PotentialEnergyDensity(J/m^3) KineticEnergyDensity(J/m^3) CrystalEnergyDensity(J/m^3)  CrystalEnergy(J)\n');
            for omegat=0:2*pi/50:16*pi
                VV=abs(dspl*dsplScale).*cos(omegat+angle(dspl)-pi/2);
                wd=0.5/Vol*transpose(VV)*DI0*VV;
                SS=abs(dspl*dsplScale).*cos(omegat+angle(dspl));
                Ekd=ow^2/2/Vol*transpose(SS)*SS;
                Vol_MoO3=4.9e-20; % Unit: m^3
                fprintf(fidout5,'%g   %g   %g   %g   %g\n',omegat/ow,wd,Ekd,wd+Ekd,(wd+Ekd)*Vol_MoO3);
            end
        end
        %==================================================
        % 与band.yaml比较
        ovlp=ph_dspl*dspl;
        ovlp=abs(ovlp);
        Epar=dot(Efield,qp(iq(i),:))/norm(qp(iq(i),:));
        Eperp=cross(Efield,qp(iq(i),:))/norm(qp(iq(i),:));
        if ow>1e11
            for j=1:nband
                if ovlp(j)>0.1 %&& minerror<1e-3   % GaP, BN这个数设为1e-3
                    fprintf(fidout2,'%g  %g  %g  %g   %g\n',norm(qp(iq(i),:)),ow/2/pi,j,ovlp(j),minerror);
                end
            end
        end
        if 1%minerror<1e-3
            fprintf(fidout,'%g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n',norm(qp(iq(i),:)),ow/2/pi,abs(real(Efield(1:3)))/norma,abs(imag(Efield(1:3)))/norma,norm(Epar),norm(Eperp),minerror);
        end
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
    end
    nwp=iwp;
    if minerror<1e-10 && (iq(i)==8)
    % 输出本征矢至nime动画
    fid2=fopen(['nime/nime',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'_',num2str(ow/2/pi/1e12,'%gTHz'),'.ascii'],'w');
    fprintf(fid2,'# Phonopy generated file for v_sim 3.6\n');
    fprintf(fid2,'%f  ',norm(lattV(1,:))*1e10);
    temp=lattV(1,:)*transpose(lattV(2,:))/norm(lattV(1,:));
    fprintf(fid2,'%f  ',temp*1e10);
    fprintf(fid2,'%f\n',sqrt(norm(lattV(2,:))^2-temp^2)*1e10);
    temp=lattV(1,:)*transpose(lattV(3,:))/norm(lattV(1,:));
    fprintf(fid2,'%f  ',temp*1e10);
    temp=cross(lattV(1,:),lattV(2,:))*transpose(lattV(3,:))/norm(cross(lattV(1,:),lattV(2,:)));
    fprintf(fid2,'%f  ',sqrt(norm(lattV(3,:))^2-norm(temp)^2)*1e10);
    fprintf(fid2,'%f\n',temp*1e10);
    for j=1:natom
        fprintf(fid2,'%f  %f  %f  B\n',atomxyz(j,1:3)*1e10);
    end
    for iwp=1:nwp
        if outputunit==1 % THz
            fprintf(fid2,'#metaData: qpt=[0.000000;0.000000;0.000000;');fprintf(fid2,'%f \\\n',wp(iwp)/2/pi/1e12);
        elseif outputunit==2 % cm-1
            fprintf(fid2,'#metaData: qpt=[0.000000;0.000000;0.000000;');fprintf(fid2,'%f \\\n',wp(iwp)/2/pi/0.03e12);
        elseif outputunit==3    % meV
            fprintf(fid2,'#metaData: qpt=[0.000000;0.000000;0.000000;');fprintf(fid2,'%f \\\n',wp(iwp)/2/pi/0.03e12*0.124);
        end
        for iatom=1:natom
            fprintf(fid2,'#');
            for idd=1:3
                fprintf(fid2,'; %f',real(dd((iatom-1)*3+idd,iwp))/sqrt(Mass(iatom)/1.66053886e-27));
            end
            for idd=1:3
                fprintf(fid2,'; %f',imag(dd((iatom-1)*3+idd,iwp))/sqrt(Mass(iatom)/1.66053886e-27));
            end
            fprintf(fid2,' \\\n');
       end
       fprintf(fid2,'# ]\n');
        % 计算振动模式的角动量Z分量
        vx=[1;0;0];
        vy=[0;1;0];
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid3,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/1e12,s,mass_s);
        elseif outputunit==2 % cm-1
            fprintf(fid3,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/0.03e12,s,mass_s);
        elseif outputunit==3    % meV
            fprintf(fid3,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/0.03e12*0.124,s,mass_s);
        end
        % 计算振动模式的角动量x分量
        vx=[0;1;0];
        vy=[0;0;1];
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid5,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/1e12,s,mass_s);
        elseif outputunit==2 % cm-1
            fprintf(fid5,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/0.03e12,s,mass_s);
        elseif outputunit==3    % meV
            fprintf(fid5,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/0.03e12*0.124,s,mass_s);
        end
        % 计算振动模式的角动量y分量
        vx=[0;0;1];
        vy=[1;0;0];
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid7,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/1e12,s,mass_s);
        elseif outputunit==2 % cm-1
            fprintf(fid7,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/0.03e12,s,mass_s);
        elseif outputunit==3    % meV
            fprintf(fid7,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g\n',iq(i),norm(qp(iq(i),:)),wp(iwp)/2/pi/0.03e12*0.124,s,mass_s);
        end
    end
     for iwp=nwp+1:nband
        fprintf(fid2,'#metaData: qpt=[0.000000;0.000000;0.000000;');fprintf(fid2,'%f \\\n',1);
        for iatom=1:natom
             fprintf(fid2,'#');
             for idd=1:3
                 fprintf(fid2,'; %f',0.01);
             end
             for idd=1:3
                 fprintf(fid2,'; %f',0);
             end
             fprintf(fid2,' \\\n');
        end
        fprintf(fid2,'# ]\n');
     end
     fclose(fid2);
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


