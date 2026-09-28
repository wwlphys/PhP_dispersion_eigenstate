fid=fopen('FORCE_CONSTANTS','r');
fscanf(fid,'%d',1);fscanf(fid,'%d',1);
ssize=dim(1)*dim(2)*dim(3);
%for i=1:natom*dim(1)*dim(2)*dim(3)
%    for j=1:natom*dim(1)*dim(2)*dim(3)
FAI=zeros(ssize*natom,ssize*natom,3,3);
for i=1:natom
    for iii=1:ssize
    for j=1:natom
        for jjj=1:ssize
        fscanf(fid,'%d',1);
        fscanf(fid,'%d',1);
        for ii=1:3
            for jj=1:3
                %FAI(i,j,ii,jj)=fscanf(fid,'%g',1);
                FAI((i-1)*ssize+iii,(j-1)*ssize+jjj,ii,jj)=fscanf(fid,'%g',1);
            end
        end
        end
    end
    end
end
fclose(fid);