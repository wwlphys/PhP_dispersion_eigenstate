function [x, y, z]=plot_cone(vstart,vend,radius)
% 画一个圆锥，起点vstart(1:3)，终点vend(1:3)，底半径radiuscone
height=norm(vend-vstart);
[x,y,z]=cylinder(radius);
[np1,np2]=size(x);
% 从圆柱变成圆锥
for i=1:np1
    for j=1:np2
        x(i,j)=x(i,j)*(1.01-z(i,j));
        y(i,j)=y(i,j)*(1.01-z(i,j));
    end
end
z=z*height; % 调整长度
% 调整圆柱体方向，先在YZ平面内旋转，再在XY平面内旋转
theta=acos((vend(3)-vstart(3))/height);%圆锥轴与Z轴的夹角
if sqrt((vend(1)-vstart(1))^2+(vend(2)-vstart(2))^2)>0
    fai=acos((vend(1)-vstart(1))/sqrt((vend(1)-vstart(1))^2+(vend(2)-vstart(2))^2));%原子间连线与X轴的夹角
else
    fai=0;
end
if vend(2)<vstart(2)
    fai=-fai;
end
for i=1:np1
    for j=1:np2
        ri=sqrt(y(i,j)^2+z(i,j)^2);
        if ri==0
            theta0=0;
        else
            theta0=acos(z(i,j)/ri);
        end
        if y(i,j)<0
            theta0=-theta0;
        end
        theta1=theta0+theta;
        y(i,j)=ri*sin(theta1);
        z(i,j)=ri*cos(theta1);
        ri=sqrt(y(i,j)^2+x(i,j)^2);
        if ri==0
            fai0=0;
        else
        fai0=acos(x(i,j)/ri);
        if y(i,j)<0
            fai0=-fai0;
        end
            fai1=(fai0-pi/2)+fai;
        end
        x(i,j)=ri*cos(fai1);
        y(i,j)=ri*sin(fai1);
    end
end
% 平移圆锥体
x=x+vstart(1);
y=y+vstart(2);
z=z+vstart(3);
