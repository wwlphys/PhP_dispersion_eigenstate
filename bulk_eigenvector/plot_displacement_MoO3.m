% 这个程序不能独自运行，只能在findeigenvector3D_largeDI.m中调用
% 1 plot x-y plane
h=figure;
plotsize(1:3)=1;
for ix=0:3
    for ilatt=0:2
        for jlatt=0:2
            plot3(real(atomxyz(:,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt+ix*20e-10),real(atomxyz(:,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt),real(atomxyz(:,3)),'o','MarkerSize',10,'Color',[0.5 0.5 0.5])
            hold on
            for iatom=1:natom
                for idd=1:3
                    position(iatom,idd)=atomxyz(iatom,idd)+dd((iatom-1)*3+idd,iwp)/sqrt(Mass(iatom)/1.66053886e-27)*6e-10*exp(ix/4*2*pi*1j);
                end
                if type(iatom)==1
                    color_data(:,iatom)=[0 0 1];
                elseif type(iatom)==2
                    color_data(:,iatom)=[1 0 0];
                elseif type(iatom)==3
                    color_data(:,iatom)=[0 1 0];
                end
                color_data(:,iatom)=color_data(:,iatom)+([1; 1; 1]-color_data(:,iatom))*0.5*(atomxyz(iatom,3)-min(atomxyz(:,3)))/(max(atomxyz(:,3)-min(atomxyz(:,3))));
                radius=0.5e-10;
                [x,y,z] = ellipsoid(position(iatom,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt+ix*20e-10,position(iatom,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt,position(iatom,3),radius,radius,radius);
                x=real(x)*plotsize(1);
                y=real(y)*plotsize(2);
                z=real(z)*plotsize(3);
                surf(x,y,z,'FaceColor',color_data(:,iatom), 'EdgeColor','none', 'FaceLighting','phong')
                if ilatt==1 && jlatt==1
                    hold on
                    %line([real(atomxyz(1,1)) real(atomxyz(1,2))],[real(atomxyz(1,1))+lattV(1,1) real(atomxyz(1,2))])
                    rectangle('position',[real(atomxyz(1,1))+lattV(1,1)*ilatt+lattV(2,1)*jlatt+ix*20e-10 real(atomxyz(1,2))+lattV(1,2)*ilatt+lattV(2,2)*jlatt lattV(1,1) lattV(2,2)],EdgeColor='g');
                    rectangle('position',[real(atomxyz(1,1))+lattV(1,1)*ilatt+lattV(2,1)*jlatt+ix*20e-10-lattV(1,1)/2 real(atomxyz(1,2))+lattV(1,2)*ilatt+lattV(2,2)*jlatt-lattV(2,2)/2 lattV(1,1)*2 lattV(2,2)*2],EdgeColor='g');
                end
            end
        end
    end
end
clear x y z;
for ii = 0:113
    x(ii+1)=min(atomxyz(:,1))+ii/100*(max(atomxyz(:,1))+8*10e-10-min(atomxyz(:,1)));
    y(ii+1)=sin(ii/100*2*pi)*lattV(2,2)*0.4;
    z(ii+1)=0;
end
x=real(x);y=real(y);z=real(z);
%plot3(x,y+lattV(2,2)*0.5,z,'--','Color','m')
%plot3(x,y+lattV(2,2)*1.5,z,'--','Color','m')
%plot3(x,y+lattV(2,2)*2.5,z,'--','Color','m')
hold off
axis equal
view(0,270)
camlight(0,0)
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
set(gcf,'OuterPosition',1e3*[-0.0062    0.0418    2.0624    1.2464])
axis off
saveas(gca,['displacement/displace_xy',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'_',num2str(iwp,'%04d'),'_',num2str(ow/2/pi/1e12,'%gTHz'),'.jpg'],'jpg')
close

% 2 plot x-z plane
h=figure;
plotsize(1:3)=1;
for ix=0:3
    for ilatt=0:2
        for klatt=0:2
            plot3(real(atomxyz(:,1)+lattV(1,1)*ilatt+ix*20e-10),real(atomxyz(:,2)+lattV(1,2)*ilatt),real(atomxyz(:,3)+lattV(3,3)*klatt),'o','MarkerSize',10,'Color',[0.5 0.5 0.5])
            hold on
            for iatom=1:natom
                for idd=1:3
                    position(iatom,idd)=atomxyz(iatom,idd)+dd((iatom-1)*3+idd,iwp)/sqrt(Mass(iatom)/1.66053886e-27)*6e-10*exp(ix/4*2*pi*1j);
                end
                if type(iatom)==1
                    color_data(:,iatom)=[0 0 1];
                elseif type(iatom)==2
                    color_data(:,iatom)=[1 0 0];
                elseif type(iatom)==3
                    color_data(:,iatom)=[0 1 0];
                end
                color_data(:,iatom)=color_data(:,iatom)-([1; 1; 1]-color_data(:,iatom))*0.5*(atomxyz(iatom,2)-max(atomxyz(:,2)))/(max(atomxyz(:,2)-min(atomxyz(:,2))));
                radius=0.5e-10;
                [x,y,z] = ellipsoid(position(iatom,1)+lattV(1,1)*ilatt+ix*20e-10,position(iatom,2)+lattV(1,2)*ilatt,position(iatom,3)+lattV(3,3)*klatt,radius,radius,radius);
                x=real(x)*plotsize(1);
                y=real(y)*plotsize(2);
                z=real(z)*plotsize(3);
                surf(x,y,z,'FaceColor',color_data(:,iatom), 'EdgeColor','none', 'FaceLighting','phong')
                if ilatt==1 && klatt==1
                    hold on
                    line([real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10) real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)+lattV(1,1) real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)+lattV(1,1) real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10) real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)],[0 0 0 0 0],[real(atomxyz(1,3)+lattV(3,3)*klatt) real(atomxyz(1,3)+lattV(3,3)*klatt) real(atomxyz(1,3)+lattV(3,3)*klatt)+lattV(3,3) real(atomxyz(1,3)+lattV(3,3)*klatt)+lattV(3,3) real(atomxyz(1,3)+lattV(3,3)*klatt)],'Color','g')
                    line([real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)-lattV(1,1)/2 real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)+lattV(1,1)*1.5 real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)+lattV(1,1)*1.5 real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)-lattV(1,1)/2 real(atomxyz(1,1)+lattV(1,1)*ilatt+ix*20e-10)-lattV(1,1)/2],[0 0 0 0 0],[real(atomxyz(1,3)+lattV(3,3)*klatt)-lattV(3,3)/4 real(atomxyz(1,3)+lattV(3,3)*klatt)-lattV(3,3)/4 real(atomxyz(1,3)+lattV(3,3)*klatt)+lattV(3,3)*1.25 real(atomxyz(1,3)+lattV(3,3)*klatt)+lattV(3,3)*1.25 real(atomxyz(1,3)+lattV(3,3)*klatt)-lattV(3,3)/4],'Color','g')
                end
            end
        end
    end
end
clear x y z;
for ii = 0:113
    x(ii+1)=min(atomxyz(:,1))+ii/100*(max(atomxyz(:,1))+8*10e-10-min(atomxyz(:,1)));
    y(ii+1)=sin(ii/100*2*pi)*lattV(3,3)*0.4;
    z(ii+1)=-sin((ii-25)/100*2*pi)*lattV(3,3)*0.4;
end
x=real(x);y=real(y);z=real(z);
%plot3(x,y,z+lattV(3,3)*0.5,'--','Color','m')
%plot3(x,y,z+lattV(3,3)*1.5,'--','Color','m')
%plot3(x,y,z+lattV(3,3)*2.5,'--','Color','m')
hold off
axis equal
view(0,180)
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
set(gcf,'OuterPosition',1e3*[-0.0062    0.0418    2.0624    1.2464])
camlight(0,0)
axis off
saveas(gca,['displacement/displace_xz',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'_',num2str(iwp,'%04d'),'_',num2str(ow/2/pi/1e12,'%gTHz'),'.jpg'],'jpg')
close
