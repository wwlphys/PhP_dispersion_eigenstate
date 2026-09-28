function [DI, colfactor]=makeDI_s(w2,miu0,epsilon0,qx,c,ndimM,DMa,Dmb,Zm,Vol,eps0,qz,thick,thick_ratio)
miuw2=miu0*w2;
emiuw2=epsilon0*miuw2;
eps_u=1;
w=sqrt(w2)/0.03e12; %wave number
delta2=0.082;omega2=806;gama2=69;delta3=0.663;omega3=1063;gama3=75; % 衬底SiO2的参数
eps_d=1.5+delta2*omega2^2/(omega2^2-w^2-1j*gama2*w)+delta3*omega3^2/(omega3^2-w^2-1j*gama3*w);%衬底SiO2相对介电常数
kuz= 1j*sqrt(qx^2-eps_u*w2/c^2);
kdz=-1j*sqrt(qx^2-eps_d*w2/c^2);
DI=zeros(ndimM*2+11,ndimM*2+17);

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

DI(2*ndimM+7,2*ndimM+2)=emiuw2-qx^2-qz^2;
DI(2*ndimM+7,2*ndimM+8)=miuw2;  % Eq.(12)

DI(2*ndimM+8,2*ndimM+2)=1;
DI(2*ndimM+8,2*ndimM+5)=1;
DI(2*ndimM+8,2*ndimM+14)=-1;  % Eq.(17)

DI(2*ndimM+9,2*ndimM+2)=exp(-1j*qz*thick);
DI(2*ndimM+9,2*ndimM+5)=exp(1j*qz*thick);
DI(2*ndimM+9,2*ndimM+17)=-1; % Eq.(18)

DI(2*ndimM+10,2*ndimM+2)=qz;
DI(2*ndimM+10,2*ndimM+5)=-qz;
DI(2*ndimM+10,2*ndimM+14)=-kuz;% Eq.(19)

DI(2*ndimM+11,2*ndimM+2)=qz*exp(-1j*qz*thick);
DI(2*ndimM+11,2*ndimM+5)=-qz*exp(1j*qz*thick);
DI(2*ndimM+11,2*ndimM+17)=-kdz;% Eq.(20)
% 缩掉X、Z分量
temp=[DI(:,1:2*ndimM) DI(:,2*ndimM+2) DI(:,2*ndimM+5) DI(:,2*ndimM+8) DI(:,2*ndimM+11) DI(:,2*ndimM+14) DI(:,2*ndimM+17)];
DI=[temp(1:2*ndimM,:);temp(2*ndimM+2,:);temp(2*ndimM+5,:);temp(2*ndimM+7:2*ndimM+11,:)];
%此时DI应为2*ndimM+7行，2*ndimM+6列
% 矩阵元素数量级相差太大，须平衡
colfactor(1:2*ndimM+6)=1;  % 记下列变换乘的系数，用于后面计算本征矢
skip=0;
if skip==0
DI(1:2*ndimM,:)=DI(1:2*ndimM,:)/max(max(abs(DI(1:2*ndimM,1:2*ndimM))));
DI(2*ndimM+1:2*ndimM+2,:)=DI(2*ndimM+1:2*ndimM+2,:)/max(max(abs(DI(2*ndimM+1:2*ndimM+2,1:2*ndimM))));
factor=max(max(abs(DI(1:2*ndimM,2*ndimM+1:2*ndimM+2))));
DI(:,2*ndimM+1:2*ndimM+2)=DI(:,2*ndimM+1:2*ndimM+2)/factor;
colfactor(2*ndimM+1:2*ndimM+2)=colfactor(2*ndimM+1:2*ndimM+2)/factor;
factor=max(max(abs(DI(2*ndimM+1:2*ndimM+2,2*ndimM+3:2*ndimM+4))));
DI(:,2*ndimM+3:2*ndimM+4)=DI(:,2*ndimM+3:2*ndimM+4)/factor;
colfactor(2*ndimM+3:2*ndimM+4)=colfactor(2*ndimM+3:2*ndimM+4)/factor;
DI(2*ndimM+3,:)=DI(2*ndimM+3,:)/sqrt(abs(DI(2*ndimM+3,2*ndimM+1))*abs(DI(2*ndimM+3,2*ndimM+3)));
DI(2*ndimM+4:2*ndimM+5,:)=DI(2*ndimM+4:2*ndimM+5,:)/abs(DI(2*ndimM+4,2*ndimM+1));
DI(2*ndimM+6:2*ndimM+7,:)=DI(2*ndimM+6:2*ndimM+7,:)/abs(DI(2*ndimM+6,2*ndimM+1));
factor=max(max(abs(DI(2*ndimM+4:2*ndimM+7,2*ndimM+5:2*ndimM+6))));
DI(:,2*ndimM+5:2*ndimM+6)=DI(:,2*ndimM+5:2*ndimM+6)/factor;
colfactor(2*ndimM+5:2*ndimM+6)=colfactor(2*ndimM+5:2*ndimM+6)/factor;
end