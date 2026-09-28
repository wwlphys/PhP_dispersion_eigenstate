% 把位移本征矢投影到vx-vy平面
clear type_s atom_s s test_norm0 test_norm chir chil;
type_s(1:ntype)=0;
atom_s(1:natom)=0;
s=0;
mass_s=0;
test_norm0=0;
test_norm=0;
chir(1:natom)=0;
chil(1:natom)=0;
for iatom=1:natom
    revec(1)=real(dd((iatom-1)*3+1,iwp));
    imvec(1)=imag(dd((iatom-1)*3+1,iwp));
    revec(2)=real(dd((iatom-1)*3+2,iwp));
    imvec(2)=imag(dd((iatom-1)*3+2,iwp));
    revec(3)=real(dd((iatom-1)*3+3,iwp));
    imvec(3)=imag(dd((iatom-1)*3+3,iwp));
    test_norm0=test_norm0+revec(1)^2+revec(2)^2+revec(3)^2+imvec(1)^2+imvec(2)^2+imvec(3)^2;
    % 投影
    pxre=dot(revec,vx);
    pxim=dot(imvec,vx);
    pyre=dot(revec,vy);
    pyim=dot(imvec,vy);
    pvec(2*iatom-1)=pxre+1j*pxim;    %  按PRB 100(2019)094303 eq.(1)构成一个列矢量
    pvec(2*iatom  )=pyre+1j*pyim;
end
% 把上述列矢量与右、左旋矢量点乘
%        pvec
for iatom=1:natom
    for j=1:natom*2
        temp(j)=right(iatom,j);
    end
    right_coef(iatom)=dot(pvec,temp);
    for j=1:natom*2
        temp(j)=left(iatom,j);
    end
    left_coef(iatom) =dot(pvec,temp);
    chir(iatom)=right_coef(iatom);
    chil(iatom)=left_coef(iatom);
    test_norm=test_norm+abs(chir(iatom))^2+abs(chil(iatom))^2;
    s=s+abs(chir(iatom))^2-abs(chil(iatom))^2; % Eq.(4) in PRB 100(2019)094303
    type_s(type(iatom))=type_s(type(iatom))+abs(chir(iatom))^2-abs(chil(iatom))^2;  % 同一种原子的才加在一起
    atom_s(iatom)=atom_s(iatom)+abs(chir(iatom))^2-abs(chil(iatom))^2;  % 同一个原子的才加在一起
    for j=1:natom*2
        temp(j)=mass_right(iatom,j);
    end
    mass_right_coef(iatom)=dot(pvec,temp);
    for j=1:natom*2
        temp(j)=mass_left(iatom,j);
    end
    mass_left_coef(iatom) =dot(pvec,temp);
    chir(iatom)=mass_right_coef(iatom);
    chil(iatom)=mass_left_coef(iatom);
    test_norm=test_norm+abs(chir(iatom))^2+abs(chil(iatom))^2;
    mass_s=mass_s+abs(chir(iatom))^2-abs(chil(iatom))^2; % 仿照Eq.(4) in PRB 100(2019)094303，只是多乘了原子质量。
end
