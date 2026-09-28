%iq=1;
DM1=zeros(nband,nband);
for kL2=1:dim(3)
for jL2=1:dim(2)
for iL2=1:dim(1)
    L2=(kL2-1)*dim(2)*dim(1)+(jL2-1)*dim(1)+iL2;
    for j=1:natom
        for j2=1:natom
            r0L2=transpose(lattV)*[iL2-1;jL2-1;kL2-1];
            r0L2=r0L2+transpose(lattV)*transpose((atomcoord(j2,:)-atomcoord(j,:)));
            normlattV(3)=sqrt(lattV(3,1)^2+lattV(3,2)^2+lattV(3,3)^2);
            % 距离近的原子之间才有FAI
            for ii=1:3
                if r0L2(ii)>=0.501*dim(ii)*norm(lattV(ii,:))
                    r0L2=r0L2-transpose(lattV(ii,:))*dim(ii);
                end
                if r0L2(ii)<=-0.501*dim(ii)*norm(lattV(ii,:))
                    r0L2=r0L2+transpose(lattV(ii,:))*dim(ii);
                end
            end
            r0L2_aaa=r0L2;
            r0L2_aab=r0L2;
            r0L2_aba=r0L2;
            r0L2_abb=r0L2;
            r0L2_baa=r0L2;
            r0L2_bab=r0L2;
            r0L2_bba=r0L2;
            r0L2_bbb=r0L2;
            % 如果两个原子在某一个方向上差超胞晶格常数的一半，则相位要平均
            if r0L2(1)>=0.499*dim(1)*norm(lattV(1,:)) && r0L2(1)<0.501*dim(1)*norm(lattV(1,:))
                r0L2_baa=r0L2_aaa-transpose(lattV(1,:))*dim(1);
                r0L2_bab=r0L2_aab-transpose(lattV(1,:))*dim(1);
                r0L2_bba=r0L2_aba-transpose(lattV(1,:))*dim(1);
                r0L2_bbb=r0L2_abb-transpose(lattV(1,:))*dim(1);
            end
            if -r0L2(1)>=0.499*dim(1)*norm(lattV(1,:)) && -r0L2(1)<0.501*dim(1)*norm(lattV(1,:))
                r0L2_baa=r0L2_aaa+transpose(lattV(1,:))*dim(1);
                r0L2_bab=r0L2_aab+transpose(lattV(1,:))*dim(1);
                r0L2_bba=r0L2_aba+transpose(lattV(1,:))*dim(1);
                r0L2_bbb=r0L2_abb+transpose(lattV(1,:))*dim(1);
            end
            if r0L2(2)>=0.499*dim(2)*norm(lattV(2,:)) && r0L2(2)<0.501*dim(2)*norm(lattV(2,:))
                r0L2_aba=r0L2_aaa-transpose(lattV(2,:))*dim(2);
                r0L2_abb=r0L2_aab-transpose(lattV(2,:))*dim(2);
                r0L2_bba=r0L2_baa-transpose(lattV(2,:))*dim(2);
                r0L2_bbb=r0L2_bab-transpose(lattV(2,:))*dim(2);
            end
            if -r0L2(2)>=0.499*dim(2)*norm(lattV(2,:)) && -r0L2(2)<0.501*dim(2)*norm(lattV(2,:))
                r0L2_aba=r0L2_aaa+transpose(lattV(2,:))*dim(2);
                r0L2_abb=r0L2_aab+transpose(lattV(2,:))*dim(2);
                r0L2_bba=r0L2_baa+transpose(lattV(2,:))*dim(2);
                r0L2_bbb=r0L2_bab+transpose(lattV(2,:))*dim(2);
            end
            if r0L2(3)>=0.499*dim(3)*norm(lattV(3,:)) && r0L2(3)<0.501*dim(3)*norm(lattV(3,:))
                r0L2_aab=r0L2_aaa-transpose(lattV(3,:))*dim(3);
                r0L2_abb=r0L2_aba-transpose(lattV(3,:))*dim(3);
                r0L2_bab=r0L2_baa-transpose(lattV(3,:))*dim(3);
                r0L2_bbb=r0L2_bba-transpose(lattV(3,:))*dim(3);
            end
            if -r0L2(3)>=0.499*dim(3)*norm(lattV(3,:)) && -r0L2(3)<0.501*dim(3)*norm(lattV(3,:))
                r0L2_aab=r0L2_aaa+transpose(lattV(3,:))*dim(3);
                r0L2_abb=r0L2_aba+transpose(lattV(3,:))*dim(3);
                r0L2_bab=r0L2_baa+transpose(lattV(3,:))*dim(3);
                r0L2_bbb=r0L2_bba+transpose(lattV(3,:))*dim(3);
            end
            expo=(exp(1j*qvec*r0L2_aaa)+exp(1j*qvec*r0L2_aab)+exp(1j*qvec*r0L2_aba)+exp(1j*qvec*r0L2_abb)+exp(1j*qvec*r0L2_baa)+exp(1j*qvec*r0L2_bab)+exp(1j*qvec*r0L2_bba)+exp(1j*qvec*r0L2_bbb))/8;
            factor=expo/sqrt(Mass(j)/1.66053886e-27)/sqrt(Mass(j2)/1.66053886e-27);
            for id=1:3
                for id2=1:3
                    DM1((j-1)*3+id,(j2-1)*3+id2)=DM1((j-1)*3+id,(j2-1)*3+id2)+FAI((j-1)*ssize+1,(j2-1)*ssize+L2,id,id2)*factor;
                end
            end
        end
    end
end
end
end
DM1=DM1*(1.602176634e-19/1.66053886e-27)/1.0e-20; % convert to standard unit
