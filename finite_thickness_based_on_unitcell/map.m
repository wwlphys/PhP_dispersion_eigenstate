clear
fid=fopen(['dispersion.txt'],'r');
num=1000000;
%iq(1:num)=0;
%diste(1:num)=0;
iqzr(1:num)=0;
ilam(1:num)=0;
%iw(1:num)=0;
%freq(1:num)=0;
qzr(1:num)=0;
qzi(1:num)=0;
mfrdiff(1:num)=0;
i=1;
iq(i)=fscanf(fid,'%d',1);
diste(i)=fscanf(fid,'%g',1);
iqzr(i)=fscanf(fid,'%d',1);
ilam(i)=fscanf(fid,'%d',1);
iw(i)=fscanf(fid,'%d',1);
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
    iw(i)=fscanf(fid,'%d',1);
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
freq=freq*0.124; % 频率单位从cm-1换成meV
nq=max(iq)-min(iq)+1;
nw=max(iw)-min(iw)+1;
mapi(1:nq,1:nw)=0;
x(1:nq,1:nw)=0;
y(1:nq,1:nw)=0;
for j=1:nw
    jj=min(iw)+j-1;
    w=min(freq)+(max(freq)-min(freq))/(nw-1)*(j-1);
    for i=1:nq
        ii=min(iq)+i-1;
        for k=1:n
            if iq(k)==ii
                mapi(i,j)=mapi(i,j)+exp(-(freq(k)-w)^2/2/7.5^2);
                x(i,j)=7/4/pi/2/diste(k);
            end
        end
        y(i,j)=min(freq)+(max(freq)-min(freq))/(nw-1)*(j-1);
    end
end
surf(x,y,mapi);
view(2)
shading interp
ylim([165 185])