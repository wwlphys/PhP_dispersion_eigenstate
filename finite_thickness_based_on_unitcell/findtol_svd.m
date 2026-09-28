function [mintol,maxtol]=findtol_svd(ncpu,DI,ndimM)
% 找临界的tol: s中从大到小第2*ndimM+11和第2*ndimM+12个，分别为maxtol和mintol
%DI=DI*sym('1');
s=svd(DI);
%s=double(s);
maxtol=s(2*ndimM+11);
mintol=s(2*ndimM+12);
%if isgpu==1
%    maxtol=gather(maxtol);
%    mintol=gather(mintol);
%    clear DI;
%end