clear
clc
c=3e8;
epsilon0=8.85e-12;
num=1000;
iq(1:num)=0;
diste(1:num)=0;
iqzr(1:num)=0;
ilam(1:num)=0;
iiw(1:num)=0;
eigenw(1:num)=0;
qzr(1:num)=0;
qzi(1:num)=0;
iphonon(1:num)=0;
project(1:num)=0;
Efieldr(1:num,1:4)=0;
PPr(1:num,1:4)=0;
Efieldvr(1:num,1:4)=0;
Efieldi(1:num,1:4)=0;
PPi(1:num,1:4)=0;
Efieldvi(1:num,1:4)=0;
tempangle(1:4)=0;
normE(1:num)=0;
tempangle2(1:3)=0;
minerror(1:num)=0;
fid=fopen('eigenvec.txt','r');
i=1;
iq(i)=fscanf(fid,'%d',1);
diste(i)=fscanf(fid,'%g',1);
iqzr(i)=fscanf(fid,'%d',1);
ilam(i)=fscanf(fid,'%d',1);
iiw(i)=fscanf(fid,'%d',1);
eigenw(i)=fscanf(fid,'%g',1)*2*pi*0.03e12/0.124;
qzr(i)=fscanf(fid,'%g',1);
qzi(i)=fscanf(fid,'%g',1);
iphonon(i)=fscanf(fid,'%d',1);
project(i)=fscanf(fid,'%g',1);
for j=1:4
    Efieldr(i,j)=fscanf(fid,'%g',1);
end
for j=1:4
    PPr(i,j)=fscanf(fid,'%g',1);
end
for j=1:4
    Efieldvr(i,j)=fscanf(fid,'%g',1);
end
for j=1:4
    Efieldi(i,j)=fscanf(fid,'%g',1);
end
for j=1:4
    PPi(i,j)=fscanf(fid,'%g',1);
end
for j=1:4
    Efieldvi(i,j)=fscanf(fid,'%g',1);
end
for j=1:4
    tempangle(j)=fscanf(fid,'%g',1);
end
normE(i)=fscanf(fid,'%g',1);
for j=1:3
    tempangle2(j)=fscanf(fid,'%g',1);
end
minerror(i)=fscanf(fid,'%g',1);
while ~feof(fid)
    i=i+1;
    iq(i)=fscanf(fid,'%d',1);
    diste(i)=fscanf(fid,'%g',1);
    iqzr(i)=fscanf(fid,'%d',1);
    ilam(i)=fscanf(fid,'%d',1);
    iiw(i)=fscanf(fid,'%d',1);
    eigenw(i)=fscanf(fid,'%g',1)*2*pi*0.03e12/0.124; % 标准单位的圆频率
    qzr(i)=fscanf(fid,'%g',1);
    qzi(i)=fscanf(fid,'%g',1);
    iphonon(i)=fscanf(fid,'%d',1);
    project(i)=fscanf(fid,'%g',1);
    for j=1:4
        Efieldr(i,j)=fscanf(fid,'%g',1);
    end
    for j=1:4
        PPr(i,j)=fscanf(fid,'%g',1);
    end
    for j=1:4
        Efieldvr(i,j)=fscanf(fid,'%g',1);
    end
    for j=1:4
        Efieldi(i,j)=fscanf(fid,'%g',1);
    end
    for j=1:4
        PPi(i,j)=fscanf(fid,'%g',1);
    end
    for j=1:4
        Efieldvi(i,j)=fscanf(fid,'%g',1);
    end
    for j=1:4
        tempangle(j)=fscanf(fid,'%g',1);
    end
    normE(i)=fscanf(fid,'%g',1);
    for j=1:3
        tempangle2(j)=fscanf(fid,'%g',1);
    end
    minerror(i)=fscanf(fid,'%g',1);
    if i>num
        disp('num should be > the number of lines of the data file')
        stop
    end
end
qz=qzr-1j*qzi;
Efield=Efieldr+1j*Efieldi;
PP=PPr+1j*PPi;
Efieldv=Efieldvr+1j*Efieldvi;
fclose(fid);
num=i;
%  读取完毕
d=10e-9
i=242
qx=diste(i);
eps_u=1;
w2=eigenw(i)^2;
w=eigenw(i)/2/pi/0.03e12; %wave number
delta2=0.082;omega2=806;gama2=69;delta3=0.663;omega3=1063;gama3=75; % 衬底SiO2的参数
eps_d=1.5+delta2*omega2^2/(omega2^2-w^2-1j*gama2*w)+delta3*omega3^2/(omega3^2-w^2-1j*gama3*w);%衬底SiO2相对介电常数
kuz=1j*sqrt(qx^2-eps_u*w2/c^2);
kdz=-1j*sqrt(qx^2-eps_d*w2/c^2);
nn=400;
xx(nn,nn)=0;
zz(nn,nn)=0;
Ezr(nn,nn)=0;
Dzr(nn,nn)=0;
for ii=1:nn
    for kk=1:nn
        x=ii/nn*3*pi/qx;
        %x=ii/nn*3*d*2;
        z=-d/2+(kk-nn/2)/nn*d*3*2;
        xx(ii,kk)=x;
        zz(ii,kk)=z;
        if z>0
            Ez=Efieldv(i,2)*exp(1j*qx*x+1j*kuz*z);
            Dz=-qx/kuz*eps_u*epsilon0*Efieldv(i,1)*exp(1j*qx*x+1j*kuz*z);
        elseif z>-d
            Ez=Efield(i,2)*exp(1j*qx*x+1j*qz(i)*z)+Efield(i,4)*exp(1j*qx*x-1j*qz(i)*z);
            Duz=-qx/kuz*eps_u*epsilon0*Efieldv(i,1); %上界面的Dz
            Ddz=-qx/kdz*eps_d*epsilon0*Efieldv(i,3); %下界面的Dz
            temp=exp(1j*qz(i)*d);
            if abs(temp-1)>1e-4
                Daz=(Duz*temp-Ddz)/(temp-1/temp);
                Dbz=(Duz/temp-Ddz)/(1/temp-temp);
            else
                Daz=Duz;
                Dbz=0;
            end
            Dz=Daz*exp(1j*qx*x+1j*qz(i)*z)+Dbz*exp(1j*qx*x-1j*qz(i)*z);
        else
            Ez=Efieldv(i,4)*exp(1j*qx*x+1j*kdz*(z+d));
            Dz=-qx/kdz*eps_d*epsilon0*Efieldv(i,3)*exp(1j*qx*x+1j*kdz*(z+d));
        end
        Ezr(ii,kk)=real(Ez);
        Dzr(ii,kk)=real(Dz);
    end
end
surf(xx,zz,Dzr);
view(0,90)
shading interp
line([xx(1,1) xx(nn,1)],[0 0],[1e100 1e100],'Color','black');
line([xx(1,1) xx(nn,1)],[-d -d],[1e100 1e100],'Color','black');