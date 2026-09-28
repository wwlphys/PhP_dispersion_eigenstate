function [mintol,maxtol]=findtol_svd_arb_q(polarization,DI,ndimM)
% 找临界的tol: s中从大到小第2*ndimM+5和第2*ndimM+6个，分别为maxtol和mintol
s=svd(DI);
if polarization ==1         % 二维材料p偏振，q须沿x轴
    maxtol=s(2*ndimM+11);
    mintol=s(2*ndimM+12);
elseif polarization == 2    % 二维材料s偏振，q须沿x轴
    maxtol=s(2*ndimM+5);
    mintol=s(2*ndimM+6);
elseif polarization == 3 % 块材
    maxtol=s(ndimM+2);   % ndimM=nband
    mintol=s(ndimM+3);
elseif polarization ==4  % 二维材料，q不可沿主轴方向
    maxtol=s(2*ndimM+17);
    mintol=s(2*ndimM+18);
end