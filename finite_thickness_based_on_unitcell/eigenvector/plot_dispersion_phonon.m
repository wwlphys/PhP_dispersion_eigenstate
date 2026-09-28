clear
fid9=fopen('dispersion_phonon_a.txt','r');
i=1;
q(i)=fscanf(fid9,'%g',1);
freq(i)=fscanf(fid9,'%g',1);
iph(i)=fscanf(fid9,'%d',1);
inproduct(i)=fscanf(fid9,'%g',1);
while ~feof(fid9)
    i=i+1;
    q(i)=fscanf(fid9,'%g',1);
    freq(i)=fscanf(fid9,'%g',1);
    iph(i)=fscanf(fid9,'%d',1);
    inproduct(i)=fscanf(fid9,'%g',1)^2;
end
scatter(q,freq,inproduct*300,iph,'filled');