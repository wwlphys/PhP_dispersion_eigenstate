function [DI, colfactor]=makeDI_arb_q(w2,miu0,epsilon0,qx,qy,c,ndimM,DMa,DMb,Zm,Vol,eps0,qz,thick,thick_ratio)
miuw2=miu0*w2;
emiuw2=epsilon0*miuw2;
eps_u=1;
w=sqrt(w2)/0.03e12; %wave number
delta2=0.082;omega2=806;gama2=69;delta3=0.663;omega3=1063;gama3=75; % 衬底SiO2的参数
eps_d=1.5+delta2*omega2^2/(omega2^2-w^2-1j*gama2*w)+delta3*omega3^2/(omega3^2-w^2-1j*gama3*w);%衬底SiO2相对介电常数
kuz=1j*sqrt(qx^2-eps_u*w2/c^2);
kdz=-1j*sqrt(qx^2-eps_d*w2/c^2);
DI=zeros(ndimM*2+22,ndimM*2+18);

DI(1:ndimM,1:ndimM)=DMa-w2*eye(ndimM);
DI(1:ndimM,2*ndimM+1:2*ndimM+3)=-1.602e-19*(Zm.');  % Eq.(13)

DI(ndimM+1:2*ndimM,ndimM+1:2*ndimM)=DMb-w2*eye(ndimM);
DI(ndimM+1:2*ndimM,2*ndimM+4:2*ndimM+6)=-1.602e-19*(Zm.'); % Eq.(14)

DI(2*ndimM+1:2*ndimM+3,1:ndimM)=1.602e-19/Vol*Zm*thick_ratio;
DI(2*ndimM+1:2*ndimM+3,2*ndimM+1:2*ndimM+3)=epsilon0*(eps0-eye(3))*thick_ratio;
DI(2*ndimM+1:2*ndimM+3,2*ndimM+7:2*ndimM+9)=-eye(3);  % Eq.(16);

DI(2*ndimM+4:2*ndimM+6,ndimM+1:2*ndimM)=1.602e-19/Vol*Zm*thick_ratio;
DI(2*ndimM+4:2*ndimM+6,2*ndimM+4:2*ndimM+6)=epsilon0*(eps0-eye(3))*thick_ratio;
DI(2*ndimM+4:2*ndimM+6,2*ndimM+10:2*ndimM+12)=-eye(3);  % Eq.(17)

DI(2*ndimM+7,2*ndimM+1)=emiuw2-qy^2-qz^2;
DI(2*ndimM+7,2*ndimM+2)=qx*qy;
DI(2*ndimM+7,2*ndimM+3)=qx*qz;
DI(2*ndimM+7,2*ndimM+7)=miuw2;  % Eq.(11)

DI(2*ndimM+8,2*ndimM+1)=qx*qy;
DI(2*ndimM+8,2*ndimM+2)=emiuw2-qx^2-qz^2;
DI(2*ndimM+8,2*ndimM+3)=qy*qz;
DI(2*ndimM+8,2*ndimM+8)=miuw2;  % Eq.(12)

DI(2*ndimM+9,2*ndimM+1)=qx*qz;
DI(2*ndimM+9,2*ndimM+2)=qy*qz;
DI(2*ndimM+9,2*ndimM+3)=emiuw2-qx^2-qy^2;
DI(2*ndimM+9,2*ndimM+9)=miuw2;  % Eq.(13)

DI(2*ndimM+10,2*ndimM+4)=emiuw2-qy^2-qz^2;
DI(2*ndimM+10,2*ndimM+5)=qx*qy;
DI(2*ndimM+10,2*ndimM+6)=-qx*qz;
DI(2*ndimM+10,2*ndimM+10)=miuw2;  % Eq.(14)

DI(2*ndimM+11,2*ndimM+4)=qx*qy;
DI(2*ndimM+11,2*ndimM+5)=emiuw2-qx^2-qz^2;
DI(2*ndimM+11,2*ndimM+6)=-qy*qz;
DI(2*ndimM+11,2*ndimM+11)=miuw2;  % Eq.(15)

DI(2*ndimM+12,2*ndimM+4)=-qx*qz;
DI(2*ndimM+12,2*ndimM+5)=-qy*qz;
DI(2*ndimM+12,2*ndimM+6)=emiuw2-qx^2-qy^2;
DI(2*ndimM+12,2*ndimM+12)=miuw2;  % Eq.(16)

DI(2*ndimM+13,2*ndimM+13)=qx;
DI(2*ndimM+13,2*ndimM+14)=qy;
DI(2*ndimM+13,2*ndimM+15)=kuz; % Eq.(17)

DI(2*ndimM+14,2*ndimM+16)=qx;
DI(2*ndimM+14,2*ndimM+17)=qy;
DI(2*ndimM+14,2*ndimM+18)=kdz; % Eq.(18)

DI(2*ndimM+15,2*ndimM+1)=1;
DI(2*ndimM+15,2*ndimM+4)=1;
DI(2*ndimM+15,2*ndimM+13)=-1;  % Eq.(19)

DI(2*ndimM+16,2*ndimM+2)=1;
DI(2*ndimM+16,2*ndimM+5)=1;
DI(2*ndimM+16,2*ndimM+14)=-1;   % eq.(20)

DI(2*ndimM+17,2*ndimM+1)=exp(-1j*qz*thick);
DI(2*ndimM+17,2*ndimM+4)=exp(1j*qz*thick);
DI(2*ndimM+17,2*ndimM+16)=-1;   % Eq.(21)

DI(2*ndimM+18,2*ndimM+2)=exp(-1j*qz*thick);
DI(2*ndimM+18,2*ndimM+5)=exp(1j*qz*thick);
DI(2*ndimM+18,2*ndimM+17)=-1;   % Eq.(22)

DI(2*ndimM+19,2*ndimM+2)=-qz;
DI(2*ndimM+19,2*ndimM+3)=qy;
DI(2*ndimM+19,2*ndimM+5)=qz;
DI(2*ndimM+19,2*ndimM+6)=qy;
DI(2*ndimM+19,2*ndimM+14)=kuz;
DI(2*ndimM+19,2*ndimM+15)=-qy;  % Eq.(23)

DI(2*ndimM+20,2*ndimM+1)=qz;
DI(2*ndimM+20,2*ndimM+3)=-qx;
DI(2*ndimM+20,2*ndimM+4)=-qz;
DI(2*ndimM+20,2*ndimM+6)=-qx;
DI(2*ndimM+20,2*ndimM+13)=-kuz;
DI(2*ndimM+20,2*ndimM+15)=qx;  % Eq.(24)

DI(2*ndimM+21,2*ndimM+2)=-qz*exp(-1j*qz*thick);
DI(2*ndimM+21,2*ndimM+3)=qy*exp(-1j*qz*thick);
DI(2*ndimM+21,2*ndimM+5)=qz*exp(1j*qz*thick);
DI(2*ndimM+21,2*ndimM+6)=qy*exp(1j*qz*thick);
DI(2*ndimM+21,2*ndimM+17)=kdz;
DI(2*ndimM+21,2*ndimM+18)=-qy; % Eq.(25)

DI(2*ndimM+22,2*ndimM+1)=qz*exp(-1j*qz*thick);
DI(2*ndimM+22,2*ndimM+3)=-qx*exp(-1j*qz*thick);
DI(2*ndimM+22,2*ndimM+4)=-qz*exp(1j*qz*thick);
DI(2*ndimM+22,2*ndimM+6)=-qx*exp(1j*qz*thick);
DI(2*ndimM+22,2*ndimM+16)=-kdz;
DI(2*ndimM+22,2*ndimM+18)=qx; % Eq.(26)
%此时DI应为2*ndimM+22行，2*ndimM+18列
% 矩阵元素数量级相差太大，须平衡
colfactor(1:2*ndimM+18)=1;  % 记下列变换乘的系数，用于后面计算本征矢
    skip=0;
if skip==0
DI(1:2*ndimM,:)=DI(1:2*ndimM,:)/max(max(abs(DI(1:2*ndimM,1:2*ndimM)))); %1
DI(2*ndimM+1:2*ndimM+6,:)=DI(2*ndimM+1:2*ndimM+6,:)/max(max(abs(DI(2*ndimM+1:2*ndimM+6,1:2*ndimM)))); %2
factor=max(max(abs(DI(1:2*ndimM,2*ndimM+1:2*ndimM+6))));
DI(:,2*ndimM+1:2*ndimM+6)=DI(:,2*ndimM+1:2*ndimM+6)/factor; %3
colfactor(2*ndimM+1:2*ndimM+6)=colfactor(2*ndimM+1:2*ndimM+6)/factor;
factor=max(max(abs(DI(2*ndimM+1:2*ndimM+6,2*ndimM+7:2*ndimM+12))));
DI(:,2*ndimM+7:2*ndimM+12)=DI(:,2*ndimM+7:2*ndimM+12)/factor; %4
colfactor(2*ndimM+7:2*ndimM+12)=colfactor(2*ndimM+7:2*ndimM+12)/factor;
DI(2*ndimM+7:2*ndimM+12,:)=DI(2*ndimM+7:2*ndimM+12,:)/sqrt(max(max(abs(DI(2*ndimM+7:2*ndimM+12,2*ndimM+1:2*ndimM+6))*max(max(abs(DI(2*ndimM+7:2*ndimM+12,2*ndimM+7:2*ndimM+12))))))); %5
%DI(2*ndimM+15:2*ndimM+18,:)=DI(2*ndimM+15:2*ndimM+18,:)/abs(DI(2*ndimM+15,2*ndimM+1)); %6
DI(2*ndimM+19:2*ndimM+22,:)=DI(2*ndimM+19:2*ndimM+22,:)/max(max(abs(DI(2*ndimM+19:2*ndimM+22,2*ndimM+1:2*ndimM+6)))); %7
factor=max(max(abs(DI(2*ndimM+19:2*ndimM+22,2*ndimM+13:2*ndimM+18))));
DI(:,2*ndimM+13:2*ndimM+18)=DI(:,2*ndimM+13:2*ndimM+18)/factor; %8
colfactor(2*ndimM+13:2*ndimM+18)=colfactor(2*ndimM+13:2*ndimM+18)/factor;
DI(2*ndimM+13:2*ndimM+14,:)=DI(2*ndimM+13:2*ndimM+14,:)/max(max(abs(DI(2*ndimM+13:2*ndimM+14,2*ndimM+13:2*ndimM+18)))); %9
DI(2*ndimM+15:2*ndimM+18,:)=DI(2*ndimM+15:2*ndimM+18,:)/sqrt(max(max(abs(DI(2*ndimM+15:2*ndimM+18,2*ndimM+1:2*ndimM+6))*max(max(abs(DI(2*ndimM+15:2*ndimM+18,2*ndimM+13:2*ndimM+18))))))); %10
end