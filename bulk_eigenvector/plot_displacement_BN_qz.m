% 这个程序不能独自运行，只能在LargeDet_minw_eigenvector_ps.m中调用
% 1 plot x-y plane
h=figure;
plotsize(1:3)=1;
for iz=0:3
    for ilatt=0:2
        for jlatt=0:2
            plot3(real(atomxyz(:,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt)+iz*12e-10,real(atomxyz(:,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt),real(atomxyz(:,3)),'o','MarkerSize',15,'Color',[0.5 0.5 0.5])
            hold on
            for iatom=1:natom
                for idd=1:3
                    position(iatom,idd)=atomxyz(iatom,idd)+dd((iatom-1)*3+idd,iwp)/sqrt(Mass(iatom)/1.66053886e-27)*6e-10/3.3*exp(iz/8*2*pi*1j);
                end
                if type(iatom)==1
                    color_data(:,iatom)=[255/255 182/255 25/255];
                elseif type(iatom)==2
                    color_data(:,iatom)=[48/255 130/255 255/255];
                elseif type(iatom)==3
                    color_data(:,iatom)=[0 1 0];
                end
                color_data(:,iatom)=color_data(:,iatom)+([1; 1; 1]-color_data(:,iatom))*0.5*(atomxyz(iatom,3)-min(atomxyz(:,3)))/(max(atomxyz(:,3)-min(atomxyz(:,3))));
                radius=1e-10/3.3;
                [x,y,z] = ellipsoid(position(iatom,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt+iz*12e-10,position(iatom,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt,position(iatom,3),radius,radius,radius);
                x=real(x)*plotsize(1);
                y=real(y)*plotsize(2);
                z=real(z)*plotsize(3);
                surf(x,y,z,'FaceColor',color_data(:,iatom), 'EdgeColor','none', 'FaceLighting','phong')
            end
        end
    end
end
clear x y z;
for ii = 0:113
    z(ii+1)=min(atomxyz(:,3))+ii/100*(max(atomxyz(:,3))+9*6e-10-min(atomxyz(:,3)));
    y(ii+1)=sin(ii/100*2*pi)*lattV(1,2)*0.4;
    x(ii+1)=0;
end
x=real(x);y=real(y);z=real(z);
plot3(x,y+lattV(1,2)*0.5,z,'--','Color','m')
plot3(x,y+lattV(1,2)*1.5,z,'--','Color','m')
plot3(x,y+lattV(1,2)*2.5,z,'--','Color','m')
hold off
axis equal
view(0,270)
camlight(0,0)
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
set(gcf,'OuterPosition',1e3*[-0.0062    0.0418    2.0624    1.2464])
saveas(gca,['displacement/displace_xy',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'_',num2str(iwp,'%04d'),'_',num2str(ow/2/pi/1e12,'%gTHz'),'.jpg'],'jpg')
close

% 2 plot x-z plane
h=figure;
plotsize(1:3)=1;
for iz=0:3
    for ilatt=0:2
        for jlatt=0:1
        for klatt=0:0
            plot3(real(atomxyz(:,3)+lattV(3,3)*klatt+iz*20e-10),real(atomxyz(:,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt),real(atomxyz(:,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt),'o','MarkerSize',10,'Color',[0.5 0.5 0.5])
            hold on
            for iatom=1:natom
                for idd=1:3
                    position(iatom,idd)=atomxyz(iatom,idd)+dd((iatom-1)*3+idd,iwp)/sqrt(Mass(iatom)/1.66053886e-27)*6e-10/3.3*exp(iz/4*2*pi*1j);
                end
                if type(iatom)==1
                    color_data(:,iatom)=[0 1 0];
                elseif type(iatom)==2
                    color_data(:,iatom)=[1 0 1];
                elseif type(iatom)==3
                    color_data(:,iatom)=[0 1 0];
                end
                color_data(:,iatom)=color_data(:,iatom)+([1; 1; 1]-color_data(:,iatom))*0.5*(atomxyz(iatom,2)-min(atomxyz(:,2)))/(max(atomxyz(:,2)-min(atomxyz(:,2))));
                radius=1e-10/3.3;
                [x,y,z] = ellipsoid(position(iatom,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt,position(iatom,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt,position(iatom,3)+lattV(3,3)*klatt+iz*20e-10,radius,radius,radius);
                x=real(x)*plotsize(1);
                y=real(y)*plotsize(2);
                z=real(z)*plotsize(3);
                surf(z,y,x,'FaceColor',color_data(:,iatom), 'EdgeColor','none', 'FaceLighting','phong')
            end
        end
        end
    end
end
clear x y z;
for ii = 0:113
    z(ii+1)=min(atomxyz(:,3))+ii/100*(max(atomxyz(:,3))+55e-10-min(atomxyz(:,3)));
    y(ii+1)=sin(ii/100*2*pi)*lattV(3,3)*0.4;
    x(ii+1)=-sin((ii-25)/100*2*pi)*lattV(3,3)*0.4;
end
x=real(x);y=real(y);z=real(z);
plot3(z,y,x+lattV(1,1)*0.5,'--','Color','m')
plot3(z,y,x+lattV(1,1)*1.5,'--','Color','m')
plot3(z,y,x+lattV(1,1)*2.5,'--','Color','m')
hold off
axis equal
view(0,180)
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
set(gcf,'OuterPosition',1e3*[-0.0062    0.0418    2.0624    1.2464])
camlight(0,0)
saveas(gca,['displacement/displace_xz',num2str(i,'%04d'),'_',num2str(iq(i),'%04d'),'_',num2str(iwp,'%04d'),'_',num2str(ow/2/pi/1e12,'%gTHz'),'.jpg'],'jpg')
close
