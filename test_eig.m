for i=1:nband
    for j=1:nband
        D(i,j)=DM(2,i,j);
    end
end
w2=eig(D)
real(sqrt(w2))/(2*pi)%*sqrt(1.6e-19/1.66e-27)/1.0e-10