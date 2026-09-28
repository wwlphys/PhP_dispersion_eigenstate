fid=fopen(['dispersion.txt'],'r');
num=1000000;
iq(1:num)=0;
diste(1:num)=0;
iqzr(1:num)=0;
ilam(1:num)=0;
freq(1:num)=0;
qzr(1:num)=0;
qzi(1:num)=0;
mfrdiff(1:num)=0;
i=1;
iq(i)=fscanf(fid,'%d',1);
diste(i)=fscanf(fid,'%g',1);
iqzr(i)=fscanf(fid,'%d',1);
ilam(i)=fscanf(fid,'%d',1);
freq(i)=fscanf(fid,'%g',1);
qzr(i)=fscanf(fid,'%g',1);
qzi(i)=fscanf(fid,'%g',1);
mfrdiff(i)=fscanf(fid,'%g',1);
while ~feof(fid)
    i=i+1;
    iq(i)=fscanf(fid,'%d',1);
    diste(i)=fscanf(fid,'%g',1);
    iqzr(i)=fscanf(fid,'%d',1);
    ilam(i)=fscanf(fid,'%d',1);
    freq(i)=fscanf(fid,'%g',1);
    qzr(i)=fscanf(fid,'%g',1);
    qzi(i)=fscanf(fid,'%g',1);
    mfrdiff(i)=fscanf(fid,'%g',1);
    if i>num
        disp('num should be > the number of lines of the data file')
        stop
    end
end
n=i;
fclose(fid);
for nn=1:10
    fid=fopen(['dispnq',num2str(nn,'%04d'),'.txt'],'w');
    for i=1:n
        fprintf(fid,'%.15g  %.15g\n',diste(i)/nn,freq(i));
    end
    fclose(fid);
end