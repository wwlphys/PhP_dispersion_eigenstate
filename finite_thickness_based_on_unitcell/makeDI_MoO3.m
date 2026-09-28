function [DI]=makeDI_MoO3(w2,miu0,epsilon0,qx,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio,isgpu)
miuw2=miu0*w2;
emiuw2=epsilon0*miuw2;
eps_u=1;
w=sqrt(w2)/0.03e12; %wave number
delta2=0.082;omega2=806;gama2=69;delta3=0.663;omega3=1063;gama3=75; % 衬底SiO2的参数
eps_d=1.5+delta2*omega2^2/(omega2^2-w^2-1j*gama2*w)+delta3*omega3^2/(omega3^2-w^2-1j*gama3*w);%衬底SiO2相对介电常数
kuz=1j*sqrt(qx^2-eps_u*w2/c^2);
kdz=-1j*sqrt(qx^2-eps_d*w2/c^2);
DI=zeros(ndimM*2+16,ndimM*2+18);

DI(1:ndimM,1:ndimM)=DMa-w2*eye(ndimM);
DI(1:ndimM,2*ndimM+1:2*ndimM+3)=-1.602e-19*(Zm.');  % Eq.(13)

DI(ndimM+1:2*ndimM,ndimM+1:2*ndimM)=DMb-w2*eye(ndimM);
DI(ndimM+1:2*ndimM,2*ndimM+4:2*ndimM+6)=-1.602e-19*(Zm.'); % Eq.(14)
% 场量在z方向依赖于qz，而原子位移不依赖于qz
%for i=1:ndimM
%    iatom=fix((i-1)/3)+1;  % 在naunitcell*nlayer个原子中是第几个原子
%    ilayer=fix((iatom-1)/natom)+1;    % 第几层的unitcell
%    iiatom=iatom-(ilayer-1)*natom;   % 这个unitcell中的第几个原子
%    zj0=unitexpdcoord((ilayer-1)*natom+iiatom,3)*lattV(i,j);
%    DI(i,2*ndimM+1:2*ndimM+3)=DI(i,2*ndimM+1:2*ndimM+3)*exp(1j*qz*zj0);
%    DI(ndimM+i,2*ndimM+1:2*ndimM+3)=DI(ndimM+i,2*ndimM+1:2*ndimM+3)*exp(1j*qz*zj0);
%end

DI(2*ndimM+1:2*ndimM+3,1:ndimM)=1.602e-19/Vol*Zm*thick_ratio;
DI(2*ndimM+1:2*ndimM+3,2*ndimM+1:2*ndimM+3)=epsilon0*(eps0-eye(3))*thick_ratio;
DI(2*ndimM+1:2*ndimM+3,2*ndimM+7:2*ndimM+9)=-eye(3);  % Eq.(16);

DI(2*ndimM+4:2*ndimM+6,ndimM+1:2*ndimM)=1.602e-19/Vol*Zm*thick_ratio;
DI(2*ndimM+4:2*ndimM+6,2*ndimM+4:2*ndimM+6)=epsilon0*(eps0-eye(3))*thick_ratio;
DI(2*ndimM+4:2*ndimM+6,2*ndimM+10:2*ndimM+12)=-eye(3);  % Eq.(17)

DI(2*ndimM+7,2*ndimM+1)=emiuw2-qz^2;
DI(2*ndimM+7,2*ndimM+3)=qx*qz;
DI(2*ndimM+7,2*ndimM+7)=miuw2;  % Eq.(11)

DI(2*ndimM+8,2*ndimM+1)=qx*qz;
DI(2*ndimM+8,2*ndimM+3)=emiuw2-qx^2;
DI(2*ndimM+8,2*ndimM+9)=miuw2;  % Eq.(12)

DI(2*ndimM+9,2*ndimM+4)=emiuw2-qz^2;
DI(2*ndimM+9,2*ndimM+6)=-qx*qz;
DI(2*ndimM+9,2*ndimM+10)=miuw2;  % Eq.(13)

DI(2*ndimM+10,2*ndimM+4)=-qx*qz;
DI(2*ndimM+10,2*ndimM+6)=emiuw2-qx^2;
DI(2*ndimM+10,2*ndimM+12)=miuw2;  % Eq.(14)

DI(2*ndimM+11,2*ndimM+13)=qx;
DI(2*ndimM+11,2*ndimM+15)=kuz; % Eq.(15)

DI(2*ndimM+12,2*ndimM+16)=qx;
DI(2*ndimM+12,2*ndimM+18)=kdz; % Eq.(16)

DI(2*ndimM+13,2*ndimM+1)=1;
DI(2*ndimM+13,2*ndimM+4)=1;
DI(2*ndimM+13,2*ndimM+13)=-1;  % Eq.(17)

DI(2*ndimM+14,2*ndimM+1)=exp(-1j*qz*thick);
DI(2*ndimM+14,2*ndimM+4)=exp(1j*qz*thick);
DI(2*ndimM+14,2*ndimM+16)=-1; % Eq.(18)

DI(2*ndimM+15,2*ndimM+1)=qz;
DI(2*ndimM+15,2*ndimM+3)=-qx;
DI(2*ndimM+15,2*ndimM+4)=-qz;
DI(2*ndimM+15,2*ndimM+6)=-qx;
DI(2*ndimM+15,2*ndimM+13)=-kuz;
DI(2*ndimM+15,2*ndimM+15)=qx;  % Eq.(19)

DI(2*ndimM+16,2*ndimM+1)=qz*exp(-1j*qz*thick);
DI(2*ndimM+16,2*ndimM+3)=-qx*exp(-1j*qz*thick);
DI(2*ndimM+16,2*ndimM+4)=-qz*exp(1j*qz*thick);
DI(2*ndimM+16,2*ndimM+6)=-qx*exp(1j*qz*thick);
DI(2*ndimM+16,2*ndimM+16)=-kdz;
DI(2*ndimM+16,2*ndimM+18)=qx; % Eq.(20)
% 缩掉Y分量
temp=[DI(:,1:2*ndimM+1) DI(:,2*ndimM+3:2*ndimM+4) DI(:,2*ndimM+6:2*ndimM+7) DI(:,2*ndimM+9:2*ndimM+10) DI(:,2*ndimM+12:2*ndimM+13) DI(:,2*ndimM+15:2*ndimM+16) DI(:,2*ndimM+18)];
DI=[temp(1:2*ndimM+1,:);temp(2*ndimM+3:2*ndimM+4,:);temp(2*ndimM+6:2*ndimM+16,:)];
%此时DI应为2*ndimM+14行，2*ndimM+12列
% 矩阵元素数量级相差太大，须平衡
skip=0;
if skip==0
if isgpu==1
    DI=gpuArray(DI);
end
DI(1:2*ndimM,:)=DI(1:2*ndimM,:)*1e-28;
DI(2*ndimM+1:2*ndimM+4,:)=DI(2*ndimM+1:2*ndimM+4,:)*1e-22;
DI(2*ndimM+5:2*ndimM+8,:)=DI(2*ndimM+5:2*ndimM+8,:)*1e-48;
if abs(qx)>1e7
    DI(2*ndimM+5:2*ndimM+8,:)=DI(2*ndimM+5:2*ndimM+8,:)*1e7/abs(qx);
end
DI(2*ndimM+9:2*ndimM+10,:)=DI(2*ndimM+9:2*ndimM+10,:)*1e-41;
if abs(qx)>1e7
    DI(2*ndimM+9:2*ndimM+10,:)=DI(2*ndimM+9:2*ndimM+10,:)*1e7/abs(qx);
end
DI(2*ndimM+11:2*ndimM+12,:)=DI(2*ndimM+11:2*ndimM+12,:)*1e-34;
DI(2*ndimM+13:2*ndimM+14,:)=DI(2*ndimM+13:2*ndimM+14,:)*1e-41;
if abs(qx)>1e7
    DI(2*ndimM+13:2*ndimM+14,:)=DI(2*ndimM+13:2*ndimM+14,:)*1e7/abs(qx);
end
DI(:,2*ndimM+1:2*ndimM+2)=DI(:,2*ndimM+1:2*ndimM+2)*1e32;
DI(:,2*ndimM+3:2*ndimM+4)=DI(:,2*ndimM+3:2*ndimM+4)*1e34;
DI(:,2*ndimM+5:2*ndimM+8)=DI(:,2*ndimM+5:2*ndimM+8)*1e24;
if abs(qx)>1e7
    DI(:,2*ndimM+5:2*ndimM+8)=DI(:,2*ndimM+5:2*ndimM+8)*sqrt(abs(qx)/1e7);
end
DI(:,2*ndimM+9:2*ndimM+12)=DI(:,2*ndimM+9:2*ndimM+12)*1e34;
end