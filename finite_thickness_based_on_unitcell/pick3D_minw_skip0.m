fid=fopen('ncompare.in','r');
ncompare=fscanf(fid,'%d',1)
dnbound=fscanf(fid,'%g',1)
upbound=fscanf(fid,'%g',1)
cqzr=fscanf(fid,'%g',1)
cqzi=fscanf(fid,'%g',1)
fclose(fid);
for iq=1:1:200
    iq
fid2=fopen(['surface_dispersion2dq',num2str(iq,'%04d'),'.txt'],'r');
fid=fopen(['surface_dispersionq',num2str(iq,'%04d'),'.txt'],'r');
if fid>0
num=4000000;
iqzr(1:num)=0;
ilam(1:num)=0;
iw(1:num)=0;
freq(1:num)=0;
qzr(1:num)=0;
qzi(1:num)=0;
mfrdiff(1:num)=0;
qp(num,3)=0;
i=1;
fscanf(fid,'%d',1);
diste=fscanf(fid,'%g',1);
iqzr(i)=fscanf(fid,'%d',1);
ilam(i)=fscanf(fid,'%d',1);
iw(i)=fscanf(fid,'%d',1);
freq(i)=fscanf(fid,'%g',1);
qzr(i)=fscanf(fid,'%g',1);
qzi(i)=fscanf(fid,'%g',1);
mfrdiff(i)=fscanf(fid,'%g',1);
%
for ii=1:3
    tq(ii)=fscanf(fid2,'%g',1);
end
qp(i,1:3)=tq(:);
fscanf(fid2,'%g',1);
while ~feof(fid)
    i=i+1;
    fscanf(fid,'%d',1);
    fscanf(fid,'%g',1);
    iqzr(i)=fscanf(fid,'%d',1);
    ilam(i)=fscanf(fid,'%d',1);
    iw(i)=fscanf(fid,'%d',1);
    freq(i)=fscanf(fid,'%g',1);
    qzr(i)=fscanf(fid,'%g',1);
    qzi(i)=fscanf(fid,'%g',1);
    mfrdiff(i)=fscanf(fid,'%g',1);
    %
    for ii=1:3
        tq(ii)=fscanf(fid2,'%g',1);
    end
    qp(i,1:3)=tq(:);
    fscanf(fid2,'%g',1);
    if i>num
        disp('num should be > the number of lines of the data file')
        stop
    end
end
fclose(fid);
fclose(fid2);
n=i;
maxiqzr=max(iqzr);
maxilam=max(ilam);
fid3=fopen(['surface_dispersion_picked_q',num2str(iq,'%04d'),'.txt'],'w');
fid4=fopen(['surface_dispersion_2d_picked_q',num2str(iq,'%04d'),'.txt'],'w');
for i=1:n
    is_min=1;
    if iqzr(i)==1 || iqzr(i)==maxiqzr || ilam(i)==1 || ilam(i)==maxilam
        is_min=0;
    end
    if abs(qzr(i))<cqzr && abs(qzi(i))<cqzi
        is_min=0;
    end
    for j=max(1,round(i-sqrt(n)*(ncompare+1))):min(n,round(i+sqrt(n)*(ncompare+1)))
        if abs(iqzr(i)-iqzr(j))<=ncompare && abs(ilam(i)-ilam(j))<=ncompare && abs(iw(i)-iw(j))<=1 && i~=j && mfrdiff(j)<=mfrdiff(i)
            is_min=0;
        end
    end
    if is_min==1
        neibor_count=0;
        for iiqzr=iqzr(i)-ncompare:iqzr(i)+ncompare
            for iilam=ilam(i)-ncompare:ilam(i)+ncompare
                if iiqzr~=iqzr(i) || iilam~=ilam(i)
                    % 在给定iiqzr,iilam的数据中找出频率与iw(i)最接近的，看是否小于mfrdiff(i)
                    mindw=1e10;
                    neibor=0;
                   for k=max(1,round(i-sqrt(n)*(ncompare+1))):min(n,round(i+sqrt(n)*(ncompare+1)))
                        if iqzr(k)==iiqzr && ilam(k)==iilam && abs(iw(k)-iw(i))<mindw
                            mindw=abs(iw(k)-iw(i));
                            neibor=k;
                        end
                    end
                    if neibor>0
                        % neibor的最近邻是不是i?
                        is_neibor=1;
                        for k=max(1,round(neibor-sqrt(n)*(ncompare+1))):min(n,round(neibor+sqrt(n)*(ncompare+1)))
                            if iqzr(k)==iqzr(i) && ilam(k)==ilam(i)  && abs(iw(k)-iw(neibor))<mindw*0.99999
                                is_neibor=0;
                            end
                        end
                        if is_neibor==1 && mfrdiff(neibor)<=mfrdiff(i)
                            is_min=0;
                        end
                        if is_neibor==1
                            neibor_count=neibor_count+1;
                        end
                    end
                end
            end
        end
        if neibor_count<(2*ncompare+1.0)^2*(3/4)
            is_min=0;
        end
    end
    if is_min==1 && freq(i)>=dnbound && freq(i)<=upbound
        fprintf(fid3,'%d  %.15g  %d  %d  %d  %.15g  %.15g  %.15g  %.15g\n',iq,diste,iqzr(i),ilam(i),iw(i),freq(i),qzr(i),qzi(i),mfrdiff(i));
        fprintf(fid4,'%.15g  %.15g  %.15g  %.15g\n',qp(i,1),qp(i,2),qp(i,3),freq(i));
    end
end
fclose(fid);
fclose(fid2);
end
end
