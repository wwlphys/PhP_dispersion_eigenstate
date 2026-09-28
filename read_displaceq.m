fid=fopen('clean_band_yaml.txt','r');
nq=fscanf(fid,'%g',1);
npath=fscanf(fid,'%g',1);
for i=1:npath
    fscanf(fid,'%g',1);
end
for i=1:3
    for j=1:3
        vb(i,j)=fscanf(fid,'%g',1)*2*pi*1e10;
    end
end
natom=fscanf(fid,'%g',1);
nband=natom*3;
for i=1:9
    fscanf(fid,'%g',1);
end
for i=1:natom
    for j=1:4
        fscanf(fid,'%g',1);
    end
    mass_au(i)=fscanf(fid,'%g',1);
end
for i=1:0       %如果band.yalm是由DFPT算的，这个数应为0；如果是有限位移法，为9.
    fscanf(fid,'%g',1);
end
freq_nac(nq,nband)=0;
ph_dspl(nq,nband,natom*3)=0;
for iq=1:nq
    for i=1:3
        qposition(iq,i)=fscanf(fid,'%g',1);
    end
    for i=1:3
        qp(iq,i)=vb(1,i)*qposition(iq,1)+vb(2,i)*qposition(iq,2)+vb(3,i)*qposition(iq,3);
    end
    fscanf(fid,'%g',1);
    for i=1:nband
        fscanf(fid,'%g',1);freq_nac(iq,i)=fscanf(fid,'%g',1);
        for j=1:natom
            fscanf(fid,'%g',1);
            for k=1:3
                re=fscanf(fid,'%g',1);im=fscanf(fid,'%g',1);
                ph_dspl(iq,i,j*3+k-3)=re+1j*im;
            end
        end
    end
end
fclose(fid);