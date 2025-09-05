%% Create figure for orientation

clear all
clc
load tipPosCrt.mat
tip_pos_x = zuSmooth{1,2}(1,:);
tip_pos_y = zuSmooth{1,2}(2,:);
tip_pos_z = zuSmooth{1,2}(3,:);
tip_ori_q0 = zuSmooth{1,2}(4,:);
tip_ori_q1 = zuSmooth{1,2}(5,:);
tip_ori_q2 = zuSmooth{1,2}(6,:);
tip_ori_q3 = zuSmooth{1,2}(7,:);
tip_times = zuSmooth{1,1};

load('actuator_motion_data_20250715_231330.mat')
spot_pos_x = -actual_coords(196:7738,1);
spot_pos_y = actual_coords(196:7738,2);

%% 根据time_real重新采样spot数据
time_real = 0.01:0.01:length(tip_pos_y)*0.01;

% 创建spot数据的原始时间轴（假设采样率相似）
spot_time_original = linspace(time_real(1), time_real(end), length(spot_pos_x));

% 使用插值将spot数据重新采样到time_real时间轴上
spot_pos_x_resampled = interp1(spot_time_original, spot_pos_x, time_real, 'linear', 'extrap');
spot_pos_y_resampled = interp1(spot_time_original, spot_pos_y, time_real, 'linear', 'extrap');

% 确保重采样后的数据与tip数据维度一致
fprintf('原始数据维度:\n');
fprintf('tip数据长度: %d\n', length(tip_pos_y));
fprintf('spot原始长度: %d\n', length(spot_pos_x));
fprintf('重采样后spot长度: %d\n', length(spot_pos_x_resampled));

%% 设计figure-8轨迹
% 基于重采样数据的中心和范围设计figure-8
center_x = mean(spot_pos_x_resampled);
center_y = mean(spot_pos_y_resampled);
range_x = max(spot_pos_x_resampled) - min(spot_pos_x_resampled);
range_y = max(spot_pos_y_resampled) - min(spot_pos_y_resampled);

% Figure-8参数
A = range_x * 0.3;  % X方向半幅（稍小于数据范围）
B = range_y * 0.15; % Y方向半幅（更小，符合figure-8特征）

% 生成figure-8轨迹点
t_figure8 = linspace(0, 2*pi, 1000);  % 高密度采样以便准确匹配
figure8_x = center_x + A * sin(t_figure8);
figure8_y = center_y + B * sin(2 * t_figure8);

fprintf('\nFigure-8轨迹参数:\n');
fprintf('中心位置: (%.2f, %.2f)\n', center_x, center_y);
fprintf('X方向半幅: %.2f mm\n', A);
fprintf('Y方向半幅: %.2f mm\n', B);
fprintf('轨迹点数: %d\n', length(figure8_x));

%% 提取与figure-8轨迹最近的采样点
tolerance = max(A, B) * 0.01;  % 容差阈值
min_spacing = max(A, B) * 0.1; % 最小间距，确保稀疏性

[extracted_indices, extracted_x, extracted_y, distances] = extract_figure8_points(...
    spot_pos_x_resampled, spot_pos_y_resampled, figure8_x, figure8_y, tolerance, min_spacing);

% 提取对应的时间
figure8_time = time_real(extracted_indices);

% 提取对应的tip位置数据
extracted_tip_x = tip_pos_x(extracted_indices);
extracted_tip_y = tip_pos_y(extracted_indices);
extracted_tip_z = tip_pos_z(extracted_indices);
extracted_tip_q0 = tip_ori_q0(extracted_indices);
extracted_tip_q1 = tip_ori_q1(extracted_indices);
extracted_tip_q2 = tip_ori_q2(extracted_indices);
extracted_tip_q3 = tip_ori_q3(extracted_indices);

% 剔除最后一个点
extracted_tip_x = extracted_tip_x(1:end-1);
extracted_tip_y = extracted_tip_y(1:end-1);
extracted_tip_z = extracted_tip_z(1:end-1);
extracted_tip_q0 = extracted_tip_q0(1:end-1);
extracted_tip_q1 = extracted_tip_q1(1:end-1);
extracted_tip_q2 = extracted_tip_q2(1:end-1);
extracted_tip_q3 = extracted_tip_q3(1:end-1);

fprintf('\n提取结果:\n');
fprintf('匹配的采样点数: %d\n', length(extracted_indices));
fprintf('平均距离误差: %.3f mm\n', mean(distances));
fprintf('最大距离误差: %.3f mm\n', max(distances));


%% 可视化结果
figure('Position', [100, 100, 900, 900]);
axis equal;
grid on;

% 子图3: 3D Tip位置轨迹
plot3(tip_pos_y, tip_pos_z, -tip_pos_x, '-', 'Color', [0.7 0.7 0.8], 'MarkerSize', 3, 'DisplayName', 'Zigzag-Scanned Tip Pose Space');
hold on;
plot3(extracted_tip_y, extracted_tip_z, -extracted_tip_x, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Tip Trajectory Reconstructed from Laser Spot Figure-8 Path');


grid on;
axis equal;
view(45, 30);

% 子图1: X-Y平面轨迹对比（原有的spot数据图）
shift = -11.26;
shift_y = -5;
spot_pos_z_resampled = shift*ones(1,size(spot_pos_x_resampled,2));
plot3(spot_pos_y_resampled, spot_pos_x_resampled + shift_y, spot_pos_z_resampled, '-', 'Color', [0.7 0.7 0.7], 'MarkerSize', 3, 'DisplayName', 'Zigzag-Scanned Laser Spot Workspace');
hold on;
figure8_z = shift*ones(1,size(figure8_x,2));
extracted_z = shift*ones(1,size(extracted_x,1));
plot3( figure8_y, figure8_x + shift_y, figure8_z,'b-', 'LineWidth', 2, 'DisplayName', 'Laser Spot Trajectory in Figure-8 Shape');
plot3( extracted_y, extracted_x + shift_y, extracted_z, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Figure-8 Trajectory from Zigzag-Scanned Workspace');

xlabel('X [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Y [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
zlabel('Z [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
title('Mapping from Laser Spot to Tip Pose', 'Fontname', 'Times New Roman', 'fontsize', 16);

% 设置legend位置为下方，并调整字体大小
legend('Location', 'southoutside', 'FontSize', 14, 'FontName', 'Times New Roman');
axis equal;
view(45, 30);
set(gca,'linewidth',2);
ax = gca;
ax.GridLineWidth = 2;
set(gca,'FontSize',16);
grid on;
axis equal;


%% 保存数据
save('figure8_extracted_data.mat', 'tip_pos_x', 'tip_pos_y', 'tip_pos_z', 'tip_times', ...
     'spot_pos_x_resampled', 'spot_pos_y_resampled', 'time_real', ...
     'spot_pos_x', 'spot_pos_y', 'figure8_x', 'figure8_y', 'extracted_indices', ...
     'extracted_x', 'extracted_y', 'extracted_tip_x', 'extracted_tip_y', 'extracted_tip_z', ...
     'extracted_tip_q0','extracted_tip_q1','extracted_tip_q2','extracted_tip_q3','figure8_time', 'distances');

fprintf('\n提取完成！\n');
fprintf('提取的采样点数: %d\n', length(extracted_indices));
fprintf('Figure-8持续时间: %.2f秒\n', figure8_time(end) - figure8_time(1));
fprintf('Tip位置变化范围:\n');
fprintf('  X: %.2f ~ %.2f mm (变化: %.2f mm)\n', min(extracted_tip_x), max(extracted_tip_x), max(extracted_tip_x) - min(extracted_tip_x));
fprintf('  Y: %.2f ~ %.2f mm (变化: %.2f mm)\n', min(extracted_tip_y), max(extracted_tip_y), max(extracted_tip_y) - min(extracted_tip_y));
fprintf('  Z: %.2f ~ %.2f mm (变化: %.2f mm)\n', min(extracted_tip_z), max(extracted_tip_z), max(extracted_tip_z) - min(extracted_tip_z));
fprintf('数据已保存至 figure8_extracted_data.mat\n');

%% 提取与figure-8轨迹最近的采样点函数
function [extracted_indices, extracted_x, extracted_y, distances] = extract_figure8_points(...
    spot_x, spot_y, figure8_x, figure8_y, tolerance, min_spacing)
    
    extracted_indices = [];
    extracted_x = [];
    extracted_y = [];
    distances = [];
    last_extracted_idx = 0;
    
    for i = 1:length(spot_x)
        % 计算当前点到figure-8轨迹的最小距离
        dist_to_curve = sqrt((figure8_x - spot_x(i)).^2 + (figure8_y - spot_y(i)).^2);
        min_dist = min(dist_to_curve);
        
        % 检查是否满足距离容差
        if min_dist <= tolerance
            % 检查与上一个提取点的间距（确保稀疏性）
            if isempty(extracted_indices) || (i - last_extracted_idx) >= min_spacing*10
                extracted_indices = [extracted_indices; i];
                extracted_x = [extracted_x; spot_x(i)];
                extracted_y = [extracted_y; spot_y(i)];
                distances = [distances; min_dist];
                last_extracted_idx = i;
            end
        end
    end
    
    % 转换为列向量
    extracted_indices = extracted_indices(:);
    extracted_x = extracted_x(:);
    extracted_y = extracted_y(:);
    distances = distances(:);
end
