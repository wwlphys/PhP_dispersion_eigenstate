function [mintol,maxtol]=findtol_svd_ps(polarization,DI,ndimM)
% 找临界的tol: s中从大到小第2*ndimM+5和第2*ndimM+6个，分别为maxtol和mintol
s=svd(DI);
if polarization ==1
    maxtol=s(2*ndimM+11);
    mintol=s(2*ndimM+12);
elseif polarization == 2
    maxtol=s(2*ndimM+5);
    mintol=s(2*ndimM+6);
elseif polarization == 3 % 块材
    maxtol=s(ndimM+2);   % ndimM=nband
    mintol=s(ndimM+3);
end