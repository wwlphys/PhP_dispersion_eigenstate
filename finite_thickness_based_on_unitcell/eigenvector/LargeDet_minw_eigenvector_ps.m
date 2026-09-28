fidin=fopen('parameters.input','r');
thick=fscanf(fidin,'%g',1)
ncpu=fscanf(fidin,'%d',1)
%parpool(ncpu);
qxi=fscanf(fidin,'%g',1)
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
else
    disp('polarization parameter error');
    stop
end
%iqout=fscanf(fidin,'%d',1) % 输出这个q点的本征矢
if is_slab==0
    thick_ratio=1;
else
    thick_ratio=lattV(3,3)/thick;
end
fclose(fidin);
wstep=(maxw-minw)/nw;
%读dispersion.txt
num=100000;
iq(1:num)=0;
diste(1:num)=0;
iqzr(1:num)=0;
ilam(1:num)=0;
iiw(1:num)=0;
eigenw(1:num)=0;
qzr(1:num)=0;
qzi(1:num)=0;
error(1:num)=0;
fid=fopen('dispersion.txt','r');
i=1;
iq(i)=fscanf(fid,'%d',1);
diste(i)=fscanf(fid,'%g',1);
iqzr(i)=fscanf(fid,'%d',1);
ilam(i)=fscanf(fid,'%d',1);
iiw(i)=fscanf(fid,'%d',1);
eigenw(i)=fscanf(fid,'%g',1);
qzr(i)=fscanf(fid,'%g',1);
qzi(i)=fscanf(fid,'%g',1);
error(i)=fscanf(fid,'%g',1);
while ~feof(fid)
    i=i+1;
    iq(i)=fscanf(fid,'%d',1);
    diste(i)=fscanf(fid,'%g',1);
    iqzr(i)=fscanf(fid,'%d',1);
    ilam(i)=fscanf(fid,'%d',1);
    iiw(i)=fscanf(fid,'%d',1);
    eigenw(i)=fscanf(fid,'%g',1);
    qzr(i)=fscanf(fid,'%g',1);
    qzi(i)=fscanf(fid,'%g',1);
    error(i)=fscanf(fid,'%g',1);
    if i>num
        disp('num should be > the number of lines of the data file')
        stop
    end
end
fclose(fid);
num=i;
%化为标准单位的圆频率
if outputunit==1
    eigenw=eigenw*2*pi;
elseif outputunit==2
    eigenw=eigenw*2*pi*0.03e12;
elseif outputunit==3
    eigenw=eigenw*2*pi*0.03e12/0.124;
else
    disp('Error');
    outputunit
end
ndimM=nband;
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
fid=fopen('eigenvec.txt','w');
fid3=fopen('Spinz_a.txt','w');
fid4=fopen('Spinz_b.txt','w');
fid5=fopen('Spinx_a.txt','w');
fid6=fopen('Spinx_b.txt','w');
fid7=fopen('Spiny_a.txt','w');
fid8=fopen('Spiny_b.txt','w');
fid9=fopen('dispersion_phonon_a.txt','w');
outcount=0;
for i=1:num
    i
    qx=sqrt(qp(iq(i),1)^2+qp(iq(i),2)^2)+1j*qxi;
    qz=qzr(i)-1j*qzi(i);
    qvec=[qx 0 qz];
    DMq
    DMa=DM1;
    qvec=[qx 0 -qz];
    DMq;
    DMb=DM1;
    w2=eigenw(i)^2;
    if polarization ==1
        [DI,colfactor]=makeDI_p(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio); % DI是2nband+14行，2nband+12列
        nrow=2*nband+14;ncol=2*nband+12;
        nfield=2; % x,z两个分量
    else
        [DI,colfactor]=makeDI_s(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio); % DI是2nband+7行，2nband+6列
        nrow=2*nband+7;ncol=2*nband+6;
        nfield=1; % y一个分量
    end
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
    minerror=min(error);       % 输出min(error)看看，min(error)太大说明这个频率不是PhP，是纯声子。
    for ideg=1:ndeg
        if error(ideg)<minerror*2
            ideg2=ideg2+1;
            vec2(:,ideg2)=vec(:,ideg);
            error(ideg2)=error(ideg);
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
        dsplb=vec2(nband+1:2*nband,ideg);   % 这是本征矢的原子振动b部分
        % 归一化
        normb=norm(dsplb);
        dsplb=dsplb/normb;
        %==============将归一化的本征矢保存===============
        iwp=iwp+1;
        wp(iwp)=eigenw(i);  % 标准单位的圆频率
        dd(:,iwp)=dspl;  
        ddb(:,iwp)=dsplb;
        %==============原子振动归一化后的电场强度=========
        Efield(1:2*nfield,iwp)=vec2(2*nband+1:2*nband+2*nfield,ideg);
        if outcount+1==204
            Efield(1:2*nfield,iwp)
        end
%        Efield(1:2)=Efield(1:2)/norma;
%        Efield(3:4)=Efield(3:4)/normb;
%       成分a、b之间有固定的比例关系，不能分别缩放
        Efield(1:2*nfield,iwp)=Efield(1:2*nfield,iwp)/sqrt(norma^2+normb^2);
        if outcount+1==204
            Efield(1:2*nfield,iwp)
        end
        PP(1:2*nfield,iwp)=vec2(2*nband+2*nfield+1:2*nband+4*nfield,ideg);
%        PP(1:2)=PP(1:2)/norma;
%        PP(3:4)=PP(3:4)/normb;
        PP(1:2*nfield,iwp)=PP(1:2*nfield,iwp)/sqrt(norma^2+normb^2);
        Efieldv(1:2*nfield,iwp)=vec2(2*nband+4*nfield+1:2*nband+6*nfield,ideg);
        Efieldv(1:2*nfield,iwp)=Efieldv(1:2*nfield,iwp)/sqrt(norma^2+normb^2);
    end
    nwp=iwp;
    %把声子本征矢做成一个矩阵
    for ieig=1:nband
        for j=1:nband
            %phonon(ieig,j)=eigvec(iq(i),ieig,j,1)+1j*eigvec(iq(i),ieig,j,2);
        end
    end
    for ii=1:nband
        for j=1:natom
            for k=1:3
                phonon(ii,j*3-3+k)=ph_dspl(ii,j*3+k-3);
            end
        end
    end
    phonon=conj(phonon);
    % 输出dd, wp，及dd与声子的内积
    for iwp=1:nwp
        for iph=1:nband
            inproduct(iph)=phonon(iph,:)*dd(:,iwp);  % 内积
        end
        inproduct=abs(inproduct);
        % 输出每个PhP在各个声子上的投影
        for iph=1:nband        
            if inproduct(iph)>1e-4
                if outputunit==1 % Hz
                    fprintf(fid9,'\n%.15g  %.15g  %d  %.15g  %g',diste(i),wp(iwp)/2/pi,iph,inproduct(iph),minerror); % 
                elseif outputunit==2 % cm-1
                    fprintf(fid9,'\n%.15g  %.15g  %d  %.15g  %g',diste(i),wp(iwp)/2/pi/0.03e12,iph,inproduct(iph),minerror); % 
                elseif outputunit==3    % meV
                    fprintf(fid9,'\n%.15g  %.15g  %d  %.15g  %g',diste(i),wp(iwp)/2/pi/0.03e12*0.124,iph,inproduct(iph),minerror); % 
                end
            end
        end
        [project,iphonon]=max(inproduct);
        if polarization ==1
            normE=sqrt(Efield(1,iwp)^2+Efield(2,iwp)^2+Efield(3,iwp)^2+Efield(4,iwp)^2);
        elseif polarization ==2
            normE=sqrt(Efield(1,iwp)^2+Efield(2,iwp)^2);
        end
        tempangle(1:2*nfield)=angle(Efield(1:2*nfield,iwp));
        for j=1:4
            if tempangle(j)<-pi/2
                tempangle(j)=tempangle(j)+2*pi;
            end
        end
        tempangle2=angle(dd(1:3,iwp));
        tempangle2(1)=tempangle2(1)-tempangle(1);
        tempangle2(3)=tempangle2(3)-tempangle(2);
        for j=1:3
            if tempangle2(j)<0
                tempangle2(j)=tempangle2(j)+2*pi;
            elseif tempangle2(j)>pi
                tempangle2(j)=tempangle2(j)-2*pi;
            end
            if tempangle2(j)<0
                tempangle2(j)=tempangle2(j)+2*pi;
            elseif tempangle2(j)>pi
                tempangle2(j)=tempangle2(j)-2*pi;
            end
        end
        if outcount+1==204
            Efield(1:2*nfield,iwp)
        end
        if outcount>=1
            fprintf(fid,"\n");
        end
        if outputunit==1 % Hz
            fprintf(fid,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %g',iq(i),qx,iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi,qzr(i),qzi(i),iphonon,project,real(Efield(1:2*nfield,iwp)),real(PP(1:2*nfield,iwp)),real(Efieldv(1:2*nfield,iwp)),imag(Efield(1:2*nfield,iwp)),imag(PP(1:2*nfield,iwp)),imag(Efieldv(1:2*nfield,iwp)),tempangle(1:2*nfield),normE,tempangle2(1:3),error(iwp)); % 
        elseif outputunit==2 % cm-1
            fprintf(fid,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %g',iq(i),qx,iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),iphonon,project,real(Efield(1:2*nfield,iwp)),real(PP(1:2*nfield,iwp)),real(Efieldv(1:2*nfield,iwp)),imag(Efield(1:2*nfield,iwp)),imag(PP(1:2*nfield,iwp)),imag(Efieldv(1:2*nfield,iwp)),tempangle(1:2*nfield),normE,tempangle2(1:3),error(iwp)); % 
        elseif outputunit==3    % meV
            fprintf(fid,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %g',iq(i),qx,iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),iphonon,project,real(Efield(1:2*nfield,iwp)),real(PP(1:2*nfield,iwp)),real(Efieldv(1:2*nfield,iwp)),imag(Efield(1:2*nfield,iwp)),imag(PP(1:2*nfield,iwp)),imag(Efieldv(1:2*nfield,iwp)),tempangle(1:2*nfield),normE,tempangle2(1:3),error(iwp)); % 
        end
        outcount=outcount+1;
        if outcount==204
       %     stop
        end
    end
    % 输出本征矢
    fid2=fopen(['nime/anime',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'.ascii'],'w');
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
%       for ieig=1:nband
%           fprintf(fid2,'#metaData: qpt=[0.000000;0.000000;0.000000;');fprintf(fid2,'%f \\\n',freq(iq(i),ieig));
%           for iatom=1:natom
%                fprintf(fid2,'#');
%                for idd=1:3
%                    fprintf(fid2,'; %f',eigvec(iq(i),ieig,(iatom-1)*3+idd,1)/sqrt(Mass(iatom)/1.66053886e-27));
%                end
%                for idd=1:3
%                    fprintf(fid2,'; %f',eigvec(iq(i),ieig,(iatom-1)*3+idd,2)/sqrt(Mass(iatom)/1.66053886e-27));
%                end
%                fprintf(fid2,' \\\n');
%           end
%           fprintf(fid2,'# ]\n');
%       end
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
       if iq(i)==200
            plot_displacement_BN_qx % 画振动图
       %else
       end
       if i==29
%           stop
       end
        % 计算振动模式的角动量Z分量
        vx=[1;0;0];
        vy=[0;1;0];
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid3,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/1e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==2 % cm-1
            fprintf(fid3,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==3    % meV
            fprintf(fid3,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        end
        % 计算振动模式的角动量x分量
        vx=[0;1;0];
        vy=[0;0;1];
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid5,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/1e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==2 % cm-1
            fprintf(fid5,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==3    % meV
            fprintf(fid5,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        end
        % 计算振动模式的角动量y分量
        vx=[0;0;1];
        vy=[1;0;0];
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid7,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/1e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==2 % cm-1
            fprintf(fid7,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==3    % meV
            fprintf(fid7,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
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
     % 输出成分b的振动图
     fid2=fopen(['nime/bnime',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'.ascii'],'w');
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
                 fprintf(fid2,'; %f',real(ddb((iatom-1)*3+idd,iwp))/sqrt(Mass(iatom)/1.66053886e-27));
             end
             for idd=1:3
                 fprintf(fid2,'; %f',imag(ddb((iatom-1)*3+idd,iwp))/sqrt(Mass(iatom)/1.66053886e-27));
             end
             fprintf(fid2,' \\\n');
        end
        fprintf(fid2,'# ]\n');
        % 计算振动模式的角动量Z分量
        vx=[1;0;0];
        vy=[0;1;0];
        dd=ddb;
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid4,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/1e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==2 % cm-1
            fprintf(fid4,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==3    % meV
            fprintf(fid4,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        end
        % 计算振动模式的角动量x分量
        vx=[0;1;0];
        vy=[0;0;1];
        dd=ddb;
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid6,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/1e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==2 % cm-1
            fprintf(fid6,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==3    % meV
            fprintf(fid6,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        end
        % 计算振动模式的角动量y分量
        vx=[0;0;1];
        vy=[1;0;0];
        dd=ddb;
        spinz_1q;
        if outputunit==1 % THz
            fprintf(fid8,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/1e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==2 % cm-1
            fprintf(fid8,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
        elseif outputunit==3    % meV
            fprintf(fid8,'%d  %g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g  %.15g  %g  %.15g  %.15g\n',iq(i),diste(i),iqzr(i),ilam(i),iiw(i),wp(iwp)/2/pi/0.03e12*0.124,qzr(i),qzi(i),s,mass_s,minerror,type_s(1:2));
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
fclose(fid);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);
fclose(fid8);
fclose(fid9);
linux=1;
if linux==1
    disp('Exit successufully!');
else
    plot(1,1);
end
