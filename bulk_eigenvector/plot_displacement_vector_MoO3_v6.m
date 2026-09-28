% 这个程序不能独自运行，只能在findeigenvector3D_largeDI.m中调用
latt_color = [0.3 0.7 0.3];
% 定义两个相位
phases = [0, pi/2]; % 相位差为π/2

for phase_idx = 1:2
    phase = phases(phase_idx);
    
    h = figure;
    plotsize(1:3) = 1;
    
    % 只绘制一个单元
    ilatt = 1;
    jlatt = 1;
    
    % 设置图形背景为白色
    set(gcf, 'Color', 'w');
    
    % 绘制平衡位置的原子
    plot3(real(atomxyz(:,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt),...
          real(atomxyz(:,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt),...
          real(atomxyz(:,3)), 'o', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5], ...
          'MarkerFaceColor', [0.8 0.8 0.8])
    hold on
    
    % 预先计算所有位移矢量，用于确定最大位移长度
    max_displacement = 0;
    for iatom = 1:natom
        displacement = zeros(1,3);
        for idd = 1:3
            displacement(idd) = real(dd((iatom-1)*3+idd, iwp) / sqrt(Mass(iatom)/1.66053886e-27) * 6e-11 * exp(1j*phase));
        end
        displacement_length = norm(displacement);
        if displacement_length > max_displacement
            max_displacement = displacement_length;
        end
    end
    
    % 设置箭头比例因子，使最大箭头长度适中
    scale_factor = 1.0/max_displacement*1.0e-10;
    
    for iatom = 1:natom
        % 计算位移矢量（考虑相位）
        displacement = zeros(1,3);
        for idd = 1:3
            displacement(idd) = real(dd((iatom-1)*3+idd, iwp) / sqrt(Mass(iatom)/1.66053886e-27) * 6e-11 * exp(1j*phase));
            if type(iatom)==1
                displacement(idd)=displacement(idd)*10;
            end
        end
        
        % 设置原子颜色
        if type(iatom) == 1
            color_data = [0.2 0.4 0.8]; % 深蓝色
            arrow_color = [0.1 0.2 0.6]; % 更深的蓝色用于箭头
        elseif type(iatom) == 2
            color_data = [0.8 0.2 0.2]; % 深红色
            arrow_color = [0.6 0.1 0.1]; % 更深的红色用于箭头
        elseif type(iatom) == 3
            color_data = [0.2 0.6 0.2]; % 深绿色
            arrow_color = [0.1 0.4 0.1]; % 更深的绿色用于箭头
        end
        
        % 调整颜色基于z坐标
        z_ratio = (atomxyz(iatom,3) - min(atomxyz(:,3))) / (max(atomxyz(:,3)) - min(atomxyz(:,3)));
        %color_data = color_data + ([0.8 0.8 0.8] - color_data) * z_ratio;
        
        % 绘制平衡位置的原子
        radius = 0.4e-10; % 稍微减小原子半径
        [x, y, z] = ellipsoid(atomxyz(iatom,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt,...
                              atomxyz(iatom,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt,...
                              atomxyz(iatom,3), radius, radius, radius);
        x = real(x) * plotsize(1);
        y = real(y) * plotsize(2);
        z = real(z) * plotsize(3);
        surf(x, y, z, 'FaceColor', color_data, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.6)
        
        % 绘制位移箭头 - 使用圆锥头部
        if norm(displacement) > 1e-12 % 只绘制有显著位移的箭头
            % 箭头起点和终点
            start_point = [atomxyz(iatom,1)+lattV(1,1)*ilatt+lattV(2,1)*jlatt,...
                          atomxyz(iatom,2)+lattV(1,2)*ilatt+lattV(2,2)*jlatt,...
                          atomxyz(iatom,3)];
            end_point = start_point + displacement * scale_factor;
            
            % 计算箭头线长度
            arrow_length = norm(displacement * scale_factor);
            
            % 绘制箭头线
            line([start_point(1), end_point(1)],...
                 [start_point(2), end_point(2)],...
                 [start_point(3), end_point(3)],...
                 'Color', arrow_color, 'LineWidth', 2.5);
            
            % 设置圆锥大小与箭头线长度成正比
            cone_length = arrow_length * 0.25; % 圆锥长度为箭头长度的25%
            cone_radius = cone_length * 0.4; % 圆锥底部半径为圆锥长度的40%
            
            % 确保圆锥有最小尺寸
            min_cone_size = 0.08e-10;
            if cone_length < min_cone_size
                cone_length = min_cone_size;
                cone_radius = min_cone_size * 0.4;
            end
            
            % 箭头方向向量
            dir_vec = displacement / norm(displacement);
            
            % 创建圆锥
            [X, Y, Z] = cylinder([0, cone_radius], 20); % 创建圆锥面
            Z = Z * cone_length; % 设置圆锥高度
            
            % 旋转圆锥以匹配箭头方向
            % 计算旋转轴和角度
            z_axis = [0, 0, 1];
            rot_axis = cross(z_axis, dir_vec);
            % 计算点积并限制在[-1, 1]内
            dot_product = dot(z_axis, dir_vec);
            if dot_product > 1
                dot_product = 1;
            elseif dot_product < -1
                dot_product = -1;
            end
            rot_angle = acos(dot_product);
            
            % 如果旋转轴长度为零（方向相同或相反），使用默认旋转
            if norm(rot_axis) < 1e-10
                if dot(z_axis, dir_vec) < 0
                    % 方向相反，旋转180度
                    R = makehgtform('axisrotate', [1, 0, 0], pi);
                else
                    % 方向相同，不需要旋转
                    R = eye(4);
                end
            else
                rot_axis = rot_axis / norm(rot_axis);
                R = makehgtform('axisrotate', rot_axis, rot_angle);
            end
            
            % 应用旋转
            for ii = 1:size(X, 1)
                for j = 1:size(X, 2)
                    point = [X(ii,j); Y(ii,j); Z(ii,j); 1];
                    rotated_point = R * point;
                    X(ii,j) = -rotated_point(1);
                    Y(ii,j) = -rotated_point(2);
                    Z(ii,j) = -rotated_point(3);
                end
            end
            
            % 修正圆锥位置：圆锥底部应该在箭头线的末端
            % 圆锥应该从箭头线的末端开始，指向位移方向
            X = X + end_point(1) + cone_length * dir_vec(1);
            Y = Y + end_point(2) + cone_length * dir_vec(2);
            Z = Z + end_point(3) + cone_length * dir_vec(3);
            
            % 绘制圆锥
            surf(X, Y, Z, 'FaceColor', arrow_color, 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
        end
    end
    
    % 绘制晶格边界
    %rectangle('Position', [real(atomxyz(1,1))+lattV(1,1)*ilatt+lattV(2,1)*jlatt,...
    %                      real(atomxyz(1,2))+lattV(1,2)*ilatt+lattV(2,2)*jlatt,...
    %                      lattV(1,1), lattV(2,2)], 'EdgeColor', [0.3 0.7 0.3], 'LineWidth', 1.5);
    for jj=1:2
        for kk=0:1
            line([real(atomxyz(1,1))+lattV(1,1)*1+lattV(2,1)*jj+lattV(3,1)*kk,       real(atomxyz(1,1))+lattV(1,1)*2+lattV(2,1)*jj+lattV(3,1)*kk],...
                 [real(atomxyz(1,2))+lattV(1,2)*1+lattV(2,2)*jj+lattV(3,2)*kk-1e-10, real(atomxyz(1,2))+lattV(1,2)*2+lattV(2,2)*jj+lattV(3,2)*kk-1e-10],...
                 [real(atomxyz(1,3))+lattV(1,3)*1+lattV(2,3)*jj+lattV(3,3)*kk,       real(atomxyz(1,3))+lattV(1,3)*2+lattV(2,3)*jj+lattV(3,3)*kk],...
                 'Color', latt_color, 'LineWidth', 1.5);
        end
    end
    for ii=1:2
        for kk=0:1
            line([real(atomxyz(1,1))+lattV(1,1)*ii+lattV(2,1)*1+lattV(3,1)*kk,       real(atomxyz(1,1))+lattV(1,1)*ii+lattV(2,1)*2+lattV(3,1)*kk],...
                 [real(atomxyz(1,2))+lattV(1,2)*ii+lattV(2,2)*1+lattV(3,2)*kk-1e-10, real(atomxyz(1,2))+lattV(1,2)*ii+lattV(2,2)*2+lattV(3,2)*kk-1e-10],...
                 [real(atomxyz(1,3))+lattV(1,3)*ii+lattV(2,3)*1+lattV(3,3)*kk,       real(atomxyz(1,3))+lattV(1,3)*ii+lattV(2,3)*2+lattV(3,3)*kk],...
                 'Color', latt_color, 'LineWidth', 1.5);
        end
    end
    for ii=1:2
        for jj=1:2
            line([real(atomxyz(1,1))+lattV(1,1)*ii+lattV(2,1)*jj+lattV(3,1)*0,       real(atomxyz(1,1))+lattV(1,1)*ii+lattV(2,1)*jj+lattV(3,1)*1],...
                 [real(atomxyz(1,2))+lattV(1,2)*ii+lattV(2,2)*jj+lattV(3,2)*0-1e-10, real(atomxyz(1,2))+lattV(1,2)*ii+lattV(2,2)*jj+lattV(3,2)*1-1e-10],...
                 [real(atomxyz(1,3))+lattV(1,3)*ii+lattV(2,3)*jj+lattV(3,3)*0,       real(atomxyz(1,3))+lattV(1,3)*ii+lattV(2,3)*jj+lattV(3,3)*1],...
                 'Color', latt_color, 'LineWidth', 1.5);
        end
    end
    % 添加坐标系指示器（xyz箭头）
    % 获取当前坐标轴范围
    ax_limits = axis;
    
    % 设置坐标系原点的位置（在图形右下角）
    origin = [ax_limits(2)*1.5, ax_limits(3)*0.8, ax_limits(5)*0.8];
    
    % 设置箭头长度（基于图形大小的比例）
    arrow_length = (ax_limits(2)-ax_limits(1)) * 0.3;
    
    % 绘制x轴（红色）
    quiver3(origin(1), origin(2), origin(3), -arrow_length, 0, 0, ...
            'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(origin(1)-arrow_length*1.1, origin(2), origin(3), 'x', ...
         'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 绘制y轴（绿色）
    quiver3(origin(1), origin(2), origin(3), 0, arrow_length, 0, ...
            'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(origin(1), origin(2)+arrow_length*1.1, origin(3), 'y', ...
         'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 绘制z轴（蓝色）
    quiver3(origin(1), origin(2), origin(3), 0, 0, arrow_length, ...
            'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(origin(1), origin(2), origin(3)+arrow_length*1.1, 'z', ...
         'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');    
    hold off
    axis equal
    view(190, -40) % 俯视图
    camlight(0, 0)
    lighting gouraud
    set(gca, 'xtick', [], 'ytick', [], 'ztick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w')
    set(gcf, 'OuterPosition', 1e3*[-0.0062, 0.0418, 2.0624, 1.2464])
    axis off
    
    % 添加标题显示频率信息
    title(sprintf('Frequency: %.2f THz, Phase: %.1fπ', ow/2/pi/1e12, phase/pi), ...
          'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
    % 保存图像，文件名中包含相位信息
    if phase_idx == 1
        phase_str = '0';
    else
        phase_str = 'pi2';
    end
    saveas(gca, ['displacement_vector/displace', '_', num2str(i,'%04d'), '_', num2str(iq(i),'%04d'), '_', num2str(iwp,'%04d'), '_', num2str(ow/2/pi/1e12,'%gTHz_phase'), phase_str, '.jpg'], 'jpg')
    close(h)
end