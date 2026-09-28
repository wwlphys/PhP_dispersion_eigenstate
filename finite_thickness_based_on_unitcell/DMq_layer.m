DM1=zeros(nband*nunitlayer,nband*nunitlayer);    %元胞*nlayer
naunitexpd=natom*nunitlayer;
for jL2=1:dim(2)
for iL2=1:dim(1)
    L2=(jL2-1)*dim(1)+iL2;
    for j=1:naunitexpd
        for j2=1:naunitexpd
            r0L2=transpose(lattV)*[iL2-1;jL2-1;0];
            r0L2=r0L2+transpose(lattV)*transpose((unitexpdcoord(j2,:)-unitexpdcoord(j,:)));
            normlattV(3)=sqrt(lattV(3,1)^2+lattV(3,2)^2+lattV(3,3)^2);
            % 距离近的原子之间才有FAI
            for i=1:2
                if r0L2(i)>=0.501*dim(i)*norm(lattV(i,:))
                    r0L2=r0L2-transpose(lattV(i,:))*dim(i);
                end
                if r0L2(i)<=-0.501*dim(i)*norm(lattV(i,:))
                    r0L2=r0L2+transpose(lattV(i,:))*dim(i);
                end
            end
            r0L2_aaa=r0L2;
            r0L2_aba=r0L2;
            r0L2_baa=r0L2;
            r0L2_bba=r0L2;
            % 如果两个原子在某一个方向上差超胞晶格常数的一半，则相位要平均
            if r0L2(1)>=0.499*dim(1)*norm(lattV(1,:)) && r0L2(1)<0.501*dim(1)*norm(lattV(1,:))
                r0L2_baa=r0L2_aaa-transpose(lattV(1,:))*dim(1);
                r0L2_bba=r0L2_aba-transpose(lattV(1,:))*dim(1);
            end
            if -r0L2(1)>=0.499*dim(1)*norm(lattV(1,:)) && -r0L2(1)<0.501*dim(1)*norm(lattV(1,:))
                r0L2_baa=r0L2_aaa+transpose(lattV(1,:))*dim(1);
                r0L2_bba=r0L2_aba+transpose(lattV(1,:))*dim(1);
            end
            if r0L2(2)>=0.499*dim(2)*norm(lattV(2,:)) && r0L2(2)<0.501*dim(2)*norm(lattV(2,:))
                r0L2_aba=r0L2_aaa-transpose(lattV(2,:))*dim(2);
                r0L2_bba=r0L2_baa-transpose(lattV(2,:))*dim(2);
            end
            if -r0L2(2)>=0.499*dim(2)*norm(lattV(2,:)) && -r0L2(2)<0.501*dim(2)*norm(lattV(2,:))
                r0L2_aba=r0L2_aaa+transpose(lattV(2,:))*dim(2);
                r0L2_bba=r0L2_baa+transpose(lattV(2,:))*dim(2);
            end
            expo=(exp(1j*qvec*r0L2_aaa)+exp(1j*qvec*r0L2_aba)+exp(1j*qvec*r0L2_baa)+exp(1j*qvec*r0L2_bba))/4;
            factor=expo/sqrt(unitexpdMass(j)/1.66053886e-27)/sqrt(unitexpdMass(j2)/1.66053886e-27);
            isupercell =fix((j -1)/natom/dim(3))+1;
            i2supercell=fix((j2-1)/natom/dim(3))+1;
            iunitcell =fix((j -(isupercell -1)*natom*dim(3)-1)/natom)*dim(1)*dim(2)+1;  %在第isupercell个超胞中的第几个元胞
            i2unitcell=fix((j2-(i2supercell-1)*natom*dim(3)-1)/natom)*dim(1)*dim(2)+1;
            iatom =(j -(isupercell -1)*natom*dim(3)-fix((iunitcell -1)/dim(1)/dim(2))*natom);  %在第isupercell个超胞中的第iunitcell个元胞中的第几个原子
            i2atom=(j2-(i2supercell-1)*natom*dim(3)-fix((i2unitcell-1)/dim(1)/dim(2))*natom);
            for id=1:3
                for id2=1:3
                    DM1((j-1)*3+id,(j2-1)*3+id2)=DM1((j-1)*3+id,(j2-1)*3+id2)+supexpdFAI((isupercell-1)*nasuper+(iatom-1)*ssize+iunitcell,(i2supercell-1)*nasuper+(i2atom-1)*ssize+i2unitcell+L2-1,id,id2)*factor;
                end
            end
        end
    end
end
end
% 场量在z方向依赖于qz，而原子位移不依赖于qz
for j=1:naunitexpd
    for j2=1:naunitexpd
        isupercell =fix((j -1)/natom/dim(3))+1;
        iunitcell =fix((j -(isupercell -1)*natom*dim(3)-1)/natom)*dim(1)*dim(2)+1;  %在第isupercell个超胞中的第几个元胞
        iatom =(j -(isupercell -1)*natom*dim(3)-fix((iunitcell -1)/dim(1)/dim(2))*natom);  %在第isupercell个超胞中的第iunitcell个元胞中的第几个原子
        ilayer=fix((j-1)/natom)+1;    % 第几层的unitcell
        zj0=unitexpdcoord((ilayer-1)*natom+iatom,3)*lattV(i,j);
        for id=1:3
            for id2=1:3
     %           DM1((j-1)*3+id,(j2-1)*3+id2)=DM1((j-1)*3+id,(j2-1)*3+id2)*exp(1j*qz*zj0);
            end
        end
    end
end

DM1=DM1*(1.602176634e-19/1.66053886e-27)/1.0e-20; % convert to standard unit
