%clear
% 如果存在path.txt，读取其中的第一行作为文件路径
% 如果不存在，则读取clearn_qpoints.txt

% 检查当前目录下是否存在path.txt文件
if exist('path.txt', 'file') == 2
    fprintf('found path.txt\n');
    
    try
        % 读取path.txt文件的第一行
        fid = fopen('path.txt', 'r');
        if fid == -1
            error('can not open path.txt');
        end
        
        target_file_path = fgetl(fid); % 读取第一行
        fclose(fid);
        
        % 检查是否成功读取到路径
        if ~ischar(target_file_path) || isempty(target_file_path)
            error('path.txt is empty or not a path');
        end
        
        fprintf('PATH read from path.txt: %s\n', target_file_path);
        
        % 检查目标文件是否存在
        if exist(target_file_path, 'file') == 2
            fprintf('reading: %s\n', target_file_path);
            % 读取目标文件内容
            fid = fopen(target_file_path, 'r');
            if fid == -1
                error('can not open: %s', target_file_path);
            end
        else
            error('file in path.txt not exist: %s', target_file_path);
        end
        
    catch ME
        fprintf('error in reading path.txt: %s\n', ME.message);
        fprintf('try to read clean_qpoints.txt\n');
        fid=fopen('clean_qpoints.txt','r');
    end
else
    fprintf('path.txt not found, trying to read clean_qpoints.txt\n');
    fid=fopen('clean_qpoints.txt','r');
end

nq=fscanf(fid,'%d',1);
natom=fscanf(fid,'%d',1);
nband=3*natom;
for i=1:3
    for j=1:3
        vb(i,j)=fscanf(fid,'%g',1)*2*pi*1e10;
    end
end
DM=zeros(nq,nband,nband);
qposition(nq,3)=0;
qp(nq,3)=0;
freq(nq,nband)=0;
eigvec(nq,nband,3*natom,2)=0;
for iq=1:nq
    if mod(iq,10000)==0
        iq
        datetime
    end
    for i=1:3
        qposition(iq,i)=fscanf(fid,'%g',1);
    end
    for i=1:3
        qp(iq,i)=vb(1,i)*qposition(iq,1)+vb(2,i)*qposition(iq,2)+vb(3,i)*qposition(iq,3);
    end
    for i=1:nband
        for j=1:nband
            re=fscanf(fid,'%g',1); im=fscanf(fid,'%g',1);
            DM(iq,i,j)=re+1j*im;
        end
    end
    for i=1:nband
        fscanf(fid,'%g',1);
        freq(iq,i)=fscanf(fid,'%g',1);
        idd=0;
        for j=1:natom
            fscanf(fid,'%g',1);
            for k=1:3
                idd=idd+1;
                eigvec(iq,i,idd,1)=fscanf(fid,'%g',1);
                eigvec(iq,i,idd,2)=fscanf(fid,'%g',1);
            end
        end
    end
end
DM=DM*(1.602176634e-19/1.66053886e-27)/1.0e-20; % convert to standard unit
