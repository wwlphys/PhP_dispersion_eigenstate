fidin=fopen('nlayer.txt','r');
nsuplayer=fscanf(fidin,'%d',1)
fclose(fidin);
nunitlayer=nsuplayer*dim(3);
nasuper=natom*ssize;
unitexpdZm=zeros(3,natom*3*nunitlayer);
unitexpdcoord=zeros(natom*nunitlayer,3);
unitexpdMass(natom*nunitlayer)=0;
supexpdcoord=zeros(nasuper*nsuplayer,3);
supexpdMass(nasuper*nsuplayer)=0;
supexpdFAI=zeros(nasuper*nsuplayer,nasuper*nsuplayer,3,3);    %超胞*nlayer
for i=1:natom
    for j=1:3
        for k=1:3
            for ilayer=1:nunitlayer
                unitexpdZm(j,(ilayer-1)*natom*3+(i-1)*3+k)=Zm(j,(i-1)*3+k);    %元胞*nlayer
            end
        end
    end
end
for ilayer=1:nunitlayer
    for iatom=1:natom
        for id=1:2
            unitexpdcoord((ilayer-1)*natom+iatom,id)=atomcoord(iatom,id);    
        end
        unitexpdcoord((ilayer-1)*natom+iatom,3)=atomcoord(iatom,3)+(ilayer-1);
        unitexpdMass((ilayer-1)*natom+iatom)=Mass(iatom); 
    end
end
for ilayer=1:nsuplayer
    for iatom=1:nasuper
        for id=1:2
            supexpdcoord((ilayer-1)*nasuper+iatom,id)=supcoord(iatom,id);    
        end
        supexpdcoord((ilayer-1)*nasuper+iatom,3)=supcoord(iatom,3)+(ilayer-1)*dim(3);
        supexpdMass((ilayer-1)*nasuper+iatom)=supMass(iatom); 
    end
end
supexpdVol=Vol*ssize*nsuplayer;
for ilayer=1:nsuplayer
    for jlayer=1:nsuplayer
        for iatom=1:nasuper
            for jatom=1:nasuper
                if abs(supexpdcoord((ilayer-1)*nasuper+iatom,3)-supexpdcoord((jlayer-1)*nasuper+jatom,3))<0.501
                    for id=1:3
                        for jd=1:3
                            supexpdFAI((ilayer-1)*nasuper+iatom,(jlayer-1)*nasuper+jatom,id,jd)=FAI(iatom,jatom,id,jd);
                        end
                    end
                end
            end
        end
    end
end
% 把supexpdFAI在XY方向不属于第一个元胞的数据删掉就是unitexpdFAI
for iunitlayer=1:nunitlayer
    isuplayer=(iunitlayer-1)/dim(3)+1;
    for junitlayer=1:nunitlayer
        jsuplayer=(junitlayer-1)/dim(3)+1;
        for iatom=1:natom
            for jatom=1:natom
                for id=1:3
                    for jd=1:3
                        unitexpdFAI((iunitlayer-1)*natom+iatom,(junitlayer-1)*natom+jatom,id,jd)=supexpdFAI(id,jd);
                    end
                end
            end
        end
    end
end
