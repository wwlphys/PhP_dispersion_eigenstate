clear
dm_preprocess
band_yaml_preprocess
%read_dm
read_displaceq
freq=freq_nac;
for iq=1:nq
    for i=1:nband
        idd=0;
        for j=1:natom
            for k=1:3
                idd=idd+1;
                eigvec(iq,i,idd,1)=real(ph_dspl(iq,i,j*3+k-3));
                eigvec(iq,i,idd,2)=imag(ph_dspl(iq,i,j*3+k-3));
            end
        end
    end
end
%  数据读取完毕
clear s;
% 看有几种原子
itype=1;
type_mass(itype)=mass_au(1);
for i=2:natom
    newtype=1;
    for j=1:i-1
        if abs(mass_au(i)-mass_au(j))<0.5
            newtype=0;
        end
    end
    if newtype==1
            itype=itype+1;
            type_mass(itype)=mass_au(i);
    end
end
ntype=itype
type_mass
for i=1:natom
    for itype=1:ntype
        if abs(mass_au(i)-type_mass(itype))<0.5
            type(i)=itype;
        end
    end
end
type
% 定义一系列natom个右旋矢量和左旋矢量，每个矢量2*natom行
for i=1:natom
    for j=1:natom*2
        right(i,j)=0;
        left(i,j)=0;
    end
    right(i,2*i-1)= 1/sqrt(2);
    right(i,2*i  )=1j/sqrt(2);
    left(i,2*i-1) =  1/sqrt(2);
    left(i,2*i  ) =-1j/sqrt(2);
end
disp('spin along q');
fid=fopen('SpinAlong_q.txt','w');
distq=0;
for iq=1:nq
    if iq==1
        qv(1)=qp(iq+1,1);qv(2)=qp(iq+1,2);qv(3)=qp(iq+1,3);
    elseif qp(iq,1)^2+qp(iq,2)^2+qp(iq,3)^2==0
        qv(1)=qp(iq-1,1);qv(2)=qp(iq-1,2);qv(3)=qp(iq-1,3);
    end
    % 在与q垂直的平面上定义vx,vy两个基矢
    temp1=[1;0;0];
    crt1=cross(temp1,qv);
    temp2=[0;1;0];
    crt2=cross(temp2,qv);
    if norm(crt1)>norm(crt2) 
        vx=crt1/norm(crt1);
    else
        vx=crt2/norm(crt2);
    end
    vy=cross(qv,vx);
    vy=vy/norm(vy);
    spinz
    if iq>1
        distq=distq+sqrt((qp(iq,1)-qp(iq-1,1))^2+(qp(iq,2)-qp(iq-1,2))^2+(qp(iq,3)-qp(iq-1,3))^2);
    end
    for ieig=1:nband
        fprintf(fid,'%g  %g  %g',distq,freq(iq,ieig),s(iq,ieig));
        for itype=1:ntype
            fprintf(fid,'  %g',type_s(iq,ieig,itype));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
disp('spin along x');
fid=fopen('SpinAlong_x.txt','w');
distq=0;
for iq=1:nq
    vx=[0;1;0];
    vy=[0;0;1];
    spinz
    if iq>1
        distq=distq+sqrt((qp(iq,1)-qp(iq-1,1))^2+(qp(iq,2)-qp(iq-1,2))^2+(qp(iq,3)-qp(iq-1,3))^2);
    end
    for ieig=1:nband
        fprintf(fid,'%g  %g  %g',distq,freq(iq,ieig),s(iq,ieig));
        for itype=1:ntype
            fprintf(fid,'  %g',type_s(iq,ieig,itype));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
disp('spin along y');
fid=fopen('SpinAlong_y.txt','w');
distq=0;
for iq=1:nq
    vx=[0;0;1];
    vy=[1;0;0];
    spinz
    if iq>1
        distq=distq+sqrt((qp(iq,1)-qp(iq-1,1))^2+(qp(iq,2)-qp(iq-1,2))^2+(qp(iq,3)-qp(iq-1,3))^2);
    end
    for ieig=1:nband
        fprintf(fid,'%g  %g  %g',distq,freq(iq,ieig),s(iq,ieig));
        for i=1:natom
            fprintf(fid,'  %g',atom_s(iq,ieig,i));
        end
        for itype=1:ntype
            fprintf(fid,'  %g',type_s(iq,ieig,itype));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
disp('spin along z');
fid=fopen('SpinAlong_z.txt','w');
distq=0;
for iq=1:nq
    vx=[1;0;0];
    vy=[0;1;0];
    spinz
    if iq>1
        distq=distq+sqrt((qp(iq,1)-qp(iq-1,1))^2+(qp(iq,2)-qp(iq-1,2))^2+(qp(iq,3)-qp(iq-1,3))^2);
    end
    for ieig=1:nband
        fprintf(fid,'%g  %g  %g',distq,freq(iq,ieig),s(iq,ieig));
        for itype=1:ntype
            fprintf(fid,'  %g',type_s(iq,ieig,itype));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
