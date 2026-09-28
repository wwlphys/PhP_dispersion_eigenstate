clear
% 选择等势面的值，你可以根据需要修改
isovalue = 35e12; 
% 把这个频率区间之外的数据删掉
lbound = 24.58e12;
ubound = 40e12;
% 读取数据，假设你的数据文件是四列的，格式为：x, y, z, V
data = load('dispersion2d.txt'); % 请替换为你的文件名
x0 = data(:, 1);
y0 = data(:, 2);
z0 = data(:, 3);
V0 = data(:, 4);
n0 = length(x0);
filter_mask = V0 <= ubound;
x_filtered = x0(filter_mask);
y_filtered = y0(filter_mask);
z_filtered = z0(filter_mask);
V_filtered = V0(filter_mask);
filter_mask2 = V_filtered >lbound;
x = x_filtered(filter_mask2);
y = y_filtered(filter_mask2);
z = z_filtered(filter_mask2);
V = V_filtered(filter_mask2);
rxy = sqrt(x.^2+y.^2+z.^2);    % 限制波矢范围
filter_mask3 = rxy < 1.3e6;
x1 = x(filter_mask3);
y1 = y(filter_mask3);
z1 = z(filter_mask3);
V1 = V(filter_mask3);
x = x1;
y = y1;
z = z1;
V = V1;
% 创建规则的网格坐标
% 这里根据你的数据范围来生成网格点，点数可以根据需要和计算资源调整
xi = linspace(min(x), max(x), 50); 
yi = linspace(min(y), max(y), 50);
zi = linspace(min(z), max(z), 50);

% 生成三维网格
[X, Y, Z] = meshgrid(xi, yi, zi);

% 将散点数据插值到规则网格上
V_grid = griddata(x, y, z, V, X, Y, Z, 'natural'); % 'natural' 方法适用于大多数情况


% 提取等值面
faces = isosurface(X, Y, Z, V_grid, isovalue);

% 绘制等值面
figure;
p = patch('Vertices', faces.vertices, 'Faces', faces.faces, ...
          'FaceColor', 'blue', 'EdgeColor', 'none');

% 设置光照和视图，让图形更美观
view(3);
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['等势面: V = ', num2str(isovalue)]);

% 添加光照
light('Position', [1, 1, 1]);
lighting gouraud;

% 设置颜色映射
colormap(jet);
%alpha(0.5)