clc
clear all
load figure8_extracted_data.mat

% 原始散点图，标出每个点的indices
figure(1);
scatter(extracted_tip_y, extracted_tip_z, 'r');
title('Original Extracted Data Points with Indices');
xlabel('Y'); ylabel('Z');
grid on;

% 在每个点旁边标出索引
for i = 1:length(extracted_tip_y)
    text(extracted_tip_y(i), extracted_tip_z(i), sprintf(' %d', i), ...
         'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold');
end

% 调整图形显示
axis equal;
hold on;

% 可选：高亮显示起点和终点
plot(extracted_tip_y(1), extracted_tip_z(1), 'go', 'MarkerSize', 12, 'LineWidth', 3);
plot(extracted_tip_y(end), extracted_tip_z(end), 'mo', 'MarkerSize', 12, 'LineWidth', 3);

% 添加图例
legend('Data Points', 'Start Point', 'End Point', 'Location', 'best');

% 手动标注留下的indices，成为时间序列的figure 8轨迹
figure8Indices = [24,20,16,12,8,4,3,2,5,9,10,14,18,21,22,27,31,35,39,43,48,47,46,45,41,37,33,29,25,24];

% 按照figure8Indices的顺序重新排列数据点
ordered_x = extracted_tip_x(figure8Indices);
ordered_y = extracted_tip_y(figure8Indices);
ordered_z = extracted_tip_z(figure8Indices);
ordered_q0 = extracted_tip_q0(figure8Indices);
ordered_q1 = extracted_tip_q1(figure8Indices);
ordered_q2 = extracted_tip_q2(figure8Indices);
ordered_q3 = extracted_tip_q3(figure8Indices);

% 创建时间向量
n_points = length(figure8Indices);
t_ordered = linspace(0, 1, n_points); % 标准化时间 0-1

% 显示重新排序后的轨迹
figure(2);
subplot(2,2,1);
plot(ordered_y, ordered_z, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('Y Position'); ylabel('Z Position');
title('Reordered Figure-8 Trajectory');
grid on;
axis equal;
hold on;
% 标记起点和终点
plot(ordered_y(1), ordered_z(1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(ordered_y(end), ordered_z(end), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('Trajectory', 'Start', 'End', 'Location', 'best');

subplot(2,2,2);
plot(t_ordered, ordered_y, 'b-o', 'MarkerSize', 4);
xlabel('Normalized Time'); ylabel('Y Position');
title('Y-axis vs Time');
grid on;

subplot(2,2,3);
plot(t_ordered, ordered_z, 'r-o', 'MarkerSize', 4);
xlabel('Normalized Time'); ylabel('Z Position');
title('Z-axis vs Time');
grid on;

subplot(2,2,4);
% 显示原始点和重新排序后的轨迹对比
scatter(extracted_tip_y, extracted_tip_z, 'ro');
hold on;
plot(ordered_y, ordered_z, 'b-', 'LineWidth', 2);
% 用数字标出重新排序后的顺序
for i = 1:length(figure8Indices)
    text(ordered_y(i), ordered_z(i), sprintf(' %d', i), ...
         'FontSize', 8, 'Color', 'green', 'FontWeight', 'bold');
end
xlabel('Y Position'); ylabel('Z Position');
title('Original Points vs Reordered Trajectory');
legend('Original Points', 'Reordered Trajectory', 'Location', 'best');
grid on;
axis equal;

ref = [ordered_x;ordered_y;ordered_z;ordered_q0;ordered_q1;ordered_q2;ordered_q3];

% 重新采样参数
total_time = 8;           % 总时间4秒
sampling_interval = 0.01; % 采样间隔0.01秒
new_time = 0:sampling_interval:total_time; % 新的时间向量
n_new_points = length(new_time); % 新的采样点数

% 原始时间向量（对应ordered数据）
original_time = linspace(0, total_time, n_points);

% 使用样条插值进行重新采样
resampled_x = interp1(original_time, ordered_x, new_time, 'spline', 'extrap');
resampled_y = interp1(original_time, ordered_y, new_time, 'spline', 'extrap');
resampled_z = interp1(original_time, ordered_z, new_time, 'spline', 'extrap');
resampled_q0 = interp1(original_time, ordered_q0, new_time, 'spline', 'extrap');
resampled_q1 = interp1(original_time, ordered_q1, new_time, 'spline', 'extrap');
resampled_q2 = interp1(original_time, ordered_q2, new_time, 'spline', 'extrap');
resampled_q3 = interp1(original_time, ordered_q3, new_time, 'spline', 'extrap');

% 重新组织ref数据
ref_resampled = [resampled_x; resampled_y; resampled_z; resampled_q0; resampled_q1; resampled_q2; resampled_q3];

% 显示重新采样后的结果
figure(3);
subplot(2,2,1);
plot(resampled_y, resampled_z, 'r-', 'LineWidth', 1.5);
hold on;
plot(ordered_y, ordered_z, 'bo-', 'MarkerSize', 6, 'LineWidth', 2);
xlabel('Y Position'); ylabel('Z Position');
title('Original vs Resampled Trajectory');
legend('Resampled (401 points)', 'Original (30 points)', 'Location', 'best');
grid on;
axis equal;

subplot(2,2,2);
plot(new_time, resampled_y, 'r-', 'LineWidth', 1.5);
hold on;
plot(original_time, ordered_y, 'bo', 'MarkerSize', 6);
xlabel('Time (s)'); ylabel('Y Position');
title('Y-axis: Original vs Resampled');
legend('Resampled', 'Original', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(new_time, resampled_z, 'r-', 'LineWidth', 1.5);
hold on;
plot(original_time, ordered_z, 'bo', 'MarkerSize', 6);
xlabel('Time (s)'); ylabel('Z Position');
title('Z-axis: Original vs Resampled');
legend('Resampled', 'Original', 'Location', 'best');
grid on;

subplot(2,2,4);
% 显示采样密度对比
stem(original_time, ones(size(original_time)), 'b', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
stem(new_time(1:10:end), 0.5*ones(size(new_time(1:10:end))), 'r', 'MarkerSize', 4);
xlabel('Time (s)'); ylabel('Sampling Points');
title('Sampling Density Comparison');
legend('Original (30 points)', 'Resampled (every 10th point)', 'Location', 'best');
grid on;
xlim([0, 4]);

% 保存重新采样后的数据
% save('figure8_resampled_trajectory.mat', 'ref_resampled', 'new_time', 'resampled_y', 'resampled_z', 'total_time', 'sampling_interval');

% 输出重新采样统计信息
fprintf('\n=== 重新采样结果 ===\n');
fprintf('原始数据点数: %d\n', n_points);
fprintf('重新采样点数: %d\n', n_new_points);
fprintf('总时间: %.2f 秒\n', total_time);
fprintf('采样间隔: %.3f 秒\n', sampling_interval);
fprintf('采样频率: %.0f Hz\n', 1/sampling_interval);
fprintf('重新采样后的轨迹长度: %.3f\n', sum(sqrt(diff(resampled_y).^2 + diff(resampled_z).^2)));

% 更新ref为重新采样后的数据
ref = ref_resampled;
%% mapping the delay-embedded ref to the low-dimensional space
% 为了对上面轨迹的第一个点做delay-embedding，考虑到这是周期性的轨迹，将这个轨迹的尾部33个time instance数据补充到前面
% mapping to the low-dimensional coordinates via chart map
for i = 1:size(ref,2)
    lowDim(:,i) = chart(ref(:,i)');
end
lowDim = [lowDim,lowDim];
for i = 1:size(ref,2)
    refReconst(:,i) = param(lowDim(:,i)');
end
% Plot trajectories for the reconstrunction accuracy comparision at the
% control observable coordinates
outdofs = [1 2 3 4 5 6 7];
figRow = length(outdofs);  % Changed to length instead of size
customFigure('subPlot',[figRow 1]); 
% Create all subplots dynamically
xyz = 'xyz'; % Create string for axis labels
for i = 1:figRow
    subplot(figRow,1,i);
    xlabel('t','Interpreter','latex');
    if i <= 3
        ylabel(sprintf('%s-axis [mm]', xyz(i)),'Interpreter','latex');
    else
        ylabel(sprintf('q%d',i-3),'Interpreter','latex');
    end
end
% Plot data for each trajectory

for i = 1:figRow
    subplot(figRow,1,i);
    plot(new_time,ref(i,:),'Linewidth',1)
    hold on
    plot(new_time,refReconst(i,:),'k:','Linewidth',2)
end


% 保存最终的ref数据
save('figure8_final_trajectory.mat', 'ref', 'new_time');
save('figure8_lowDim_trajectory.mat', 'lowDim');

fprintf('\nmapping完成！新的ref矩阵尺寸: %dx%d\n', size(ref));
%% Create timeseries objects for Simulink From Workspace block
% Create timeseries for each quaternion component
ts_sensor2Pos1 = timeseries(resampled_x, new_time, 'Name', 'sensor2_Pos_x');
ts_sensor2Pos1.DataInfo.Units = '';
ts_sensor2Pos1.DataInfo.Interpolation = tsdata.interpolation('linear');

ts_sensor2Pos2 = timeseries(resampled_y, new_time, 'Name', 'sensor2_Pos_y');
ts_sensor2Pos2.DataInfo.Units = '';
ts_sensor2Pos2.DataInfo.Interpolation = tsdata.interpolation('linear');

ts_sensor2Pos3 = timeseries(resampled_z, new_time, 'Name', 'sensor2_Pos_z');
ts_sensor2Pos3.DataInfo.Units = '';
ts_sensor2Pos3.DataInfo.Interpolation = tsdata.interpolation('linear');

ts_sensor2Quater1 = timeseries(resampled_q0, new_time, 'Name', 'sensor2_Quaternion_w');
ts_sensor2Quater1.DataInfo.Units = '';
ts_sensor2Quater1.DataInfo.Interpolation = tsdata.interpolation('linear');

ts_sensor2Quater2 = timeseries(resampled_q1, new_time, 'Name', 'sensor2_Quaternion_x');
ts_sensor2Quater2.DataInfo.Units = '';
ts_sensor2Quater2.DataInfo.Interpolation = tsdata.interpolation('linear');

ts_sensor2Quater3 = timeseries(resampled_q2, new_time, 'Name', 'sensor2_Quaternion_y');
ts_sensor2Quater3.DataInfo.Units = '';
ts_sensor2Quater3.DataInfo.Interpolation = tsdata.interpolation('linear');

ts_sensor2Quater4 = timeseries(resampled_q3, new_time, 'Name', 'sensor2_Quaternion_z');
ts_sensor2Quater4.DataInfo.Units = '';
ts_sensor2Quater4.DataInfo.Interpolation = tsdata.interpolation('linear');

% Save timeseries to .mat file for Simulink
save('ref_lowDim_Timeseries.mat', 'ts_sensor2Pos1', 'ts_sensor2Pos2', 'ts_sensor2Pos3', 'ts_sensor2Quater1', 'ts_sensor2Quater2', ...
     'ts_sensor2Quater3', 'ts_sensor2Quater4');
