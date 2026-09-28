% 把位移本征矢投影到vx-vy平面
if iq==1
    type_s(1:nq,1:nband,1:ntype)=0;
    atom_s(1:nq,1:nband,1:natom)=0;
    s(1:nq,1:nband)=0;
    test_norm0(1:nq,1:nband)=0;
    test_norm(1:nq,1:nband)=0;
    chir(1:nq,1:nband,1:natom)=0;
    chil(1:nq,1:nband,1:natom)=0;
end
    for ieig=1:nband
        for i=1:natom
            revec(1)=eigvec(iq,ieig,(i-1)*3+1,1);
            imvec(1)=eigvec(iq,ieig,(i-1)*3+1,2);
            revec(2)=eigvec(iq,ieig,(i-1)*3+2,1);
            imvec(2)=eigvec(iq,ieig,(i-1)*3+2,2);
            revec(3)=eigvec(iq,ieig,(i-1)*3+3,1);
            imvec(3)=eigvec(iq,ieig,(i-1)*3+3,2);
            test_norm0(iq,ieig)=test_norm0(iq,ieig)+revec(1)^2+revec(2)^2+revec(3)^2+imvec(1)^2+imvec(2)^2+imvec(3)^2;
            % 投影
            pxre=dot(revec,vx);
            pxim=dot(imvec,vx);
            pyre=dot(revec,vy);
            pyim=dot(imvec,vy);
            pvec(2*i-1)=pxre+1j*pxim;    %  按PRB 100(2019)094303 eq.(1)构成一个列矢量
            pvec(2*i  )=pyre+1j*pyim;
        end
        % 把上述列矢量与右、左旋矢量点乘
%        pvec
        for i=1:natom
            for j=1:natom*2
                temp(j)=right(i,j);
            end
%            temp
            right_coef(i)=dot(pvec,temp);
            for j=1:natom*2
                temp(j)=left(i,j);
            end
%            temp
            left_coef(i) =dot(pvec,temp);
            chir(iq,ieig,i)=right_coef(i);
            chil(iq,ieig,i)=left_coef(i);
            test_norm(iq,ieig)=test_norm(iq,ieig)+abs(chir(iq,ieig,i))^2+abs(chil(iq,ieig,i))^2;
            s(iq,ieig)=s(iq,ieig)+abs(chir(iq,ieig,i))^2-abs(chil(iq,ieig,i))^2; % Eq.(4) in PRB 100(2019)094303
            type_s(iq,ieig,type(i))=type_s(iq,ieig,type(i))+abs(chir(iq,ieig,i))^2-abs(chil(iq,ieig,i))^2;  % 同一种原子的才加在一起
            atom_s(iq,ieig,i)=atom_s(iq,ieig,i)+abs(chir(iq,ieig,i))^2-abs(chil(iq,ieig,i))^2;  % 同一个原子的才加在一起
%            stop
        end
    end
