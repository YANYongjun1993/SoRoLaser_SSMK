%% load new control data
clearvars;
clc
clear all
close all
% Load new controller data with random inputs
load 202505252158.mat
load SSM_model.mat
sensor1_Quater1 = get(data,"sensor1_Quater1");
time = sensor1_Quater1.Values.Time;
Sensor1Quater1 = sensor1_Quater1.Values.data;
sensor1_Quater2 = get(data,"sensor1_Quater2");
Sensor1Quater2 = sensor1_Quater2.Values.data;
sensor1_Quater3 = get(data,"sensor1_Quater3");
Sensor1Quater3 = sensor1_Quater3.Values.data;
sensor1_Quater4 = get(data,"sensor1_Quater4");
Sensor1Quater4 = sensor1_Quater4.Values.data;
sensor1_x = get(data,"sensor1_x");
Sensor1X = sensor1_x.Values.data;
sensor1_y = get(data,"sensor1_y");
Sensor1Y = sensor1_y.Values.data;
sensor1_z = get(data,"sensor1_z");
Sensor1Z = sensor1_z.Values.data;
sensor2_Quater1 = get(data,"sensor2_Quater1");
Sensor2Quater1 = sensor2_Quater1.Values.data;
sensor2_Quater2 = get(data,"sensor2_Quater2");
Sensor2Quater2 = sensor2_Quater2.Values.data;
sensor2_Quater3 = get(data,"sensor2_Quater3");
Sensor2Quater3 = sensor2_Quater3.Values.data;
sensor2_Quater4 = get(data,"sensor2_Quater4");
Sensor2Quater4 = sensor2_Quater4.Values.data;
sensor2_x = get(data,"sensor2_x");
Sensor2X = sensor2_x.Values.data;
sensor2_y = get(data,"sensor2_y");
Sensor2Y = sensor2_y.Values.data;
sensor2_z = get(data,"sensor2_z");
Sensor2Z = sensor2_z.Values.data;
% concatenating the position and quaterion
Sensor2 = [Sensor2X';Sensor2Y';Sensor2Z';Sensor2Quater1';Sensor2Quater2';Sensor2Quater3';Sensor2Quater4'];

flag = get(data,"flag");
flagStatus = flag.Values.data;

Motor1_pos = get(data,"Motor1_pos_req_qc_mode2");
Motor1Pos = Motor1_pos.Values.data;
Motor2_pos = get(data,"Motor2_pos_req_qc_mode2");
Motor2Pos = Motor2_pos.Values.data;
Motor3_pos = get(data,"Motor3_pos_req_qc_mode2");
Motor3Pos = Motor3_pos.Values.data;
Motor4_pos = get(data,"Motor4_pos_req_qc_mode2");
Motor4Pos = Motor4_pos.Values.data;

uMotorPos = [Motor1Pos';Motor2Pos';Motor3Pos';Motor4Pos'];
%%
B = ssm_model.B;
Br = ssm_model.Br;
Rauton = ssm_model.RDInfoCtrl.reducedDynamics.map;
Monom = ssm_model.RDInfoCtrl.reducedDynamics.phi;
R = @(t,x,u) Rauton(x) + B*u + Br*kron(Monom(x),u);
%% setup matrix
zu = {};
keptIdxs = 10960:30950;
y_eq = Sensor2(:,31266);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor2(:, keptIdxs) - y_eq);
u_eq = uMotorPos(:, 31266);
uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorPos(:, keptIdxs) - u_eq; 
%% Multi-stage smoothing for reliable differential calculus
% Parameters
windowSize = 101;  % Must be odd
polynomialOrder = 2;  % Balance between smoothing and feature preservation
medianFilterSize = 7;  % For removing spikes
butterworthOrder = 4;  % Butterworth filter order
cutoffFreq = 0.01;     % Normalized cutoff frequency (0 to 1)
  
for iDim = 1:size(Sensor2,1)
    % 1. Apply median filter to remove spikes
    temp_data = medfilt1(zu{1,2} (iDim,:), medianFilterSize);
    
    % 2. Apply Butterworth filter
    [b, a] = butter(butterworthOrder, cutoffFreq, 'low');
    temp_data = filtfilt(b, a, temp_data);
    
    % 3. Final smoothing with Savitzky-Golay
    Sensor2SmoothTemp(iDim,:) = sgolayfilt(temp_data, polynomialOrder, windowSize);
end

% Store smoothed data
zuSmooth{1,2} = Sensor2SmoothTemp;
zuSmooth{1,1} = zu{1,1};

% Plot trajectories for the filter effect visualization
outdofs = [1 2 3];
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
    plot(zu{1,1},zu{1,2}(outdofs(i),:),'Linewidth',1)
    hold on
    plot(zuSmooth{1,1},zuSmooth{1,2}(outdofs(i),:),'k:','Linewidth',2)
end

Yu = zuSmooth; 
%% check if the controller trajectories are on the manifold, v then w mapping
xData = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
yDataRemapped = liftTrajectories(ssm_model.IMInfoCtrl, xData); % Project up

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
    plot(yDataRemapped{1,1},yDataRemapped{1,2}(outdofs(i),:),'Linewidth',1)
    hold on
    plot(Yu{1,1},Yu{1,2}(outdofs(i),:),'k:','Linewidth',2)
end

%% Validate position

figure('Position', [100, 100, 600, 600]);
axis equal;
grid on;
xlabel('X [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Y [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
zlabel('Z [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
title('Tip Orientation Over Time', 'Fontname', 'Times New Roman', 'fontsize', 16);
hold on;

plot3(Yu{1,2}(2,1:10:end), Yu{1,2}(3,1:10:end), -Yu{1,2}(1,1:10:end),'-','Color', [0.8 0.2 0.2], 'LineWidth', 2);
hold on;
plot3(yDataRemapped{1,2}(2,1:10:end), yDataRemapped{1,2}(3,1:10:end), -yDataRemapped{1,2}(1,1:10:end),':','Color', [0.1 0.3 0.8], 'LineWidth', 3);

zlim([-1,1]);
view([-30 40]);
set(gca,'linewidth',2);
ax = gca;
ax.GridLineWidth = 2;
set(gca,'FontSize',16);
grid on;
axis equal;

%% orientation
figure('Position', [100, 100, 600, 600]);
axis equal;
grid on;
xlabel('X [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Y [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
zlabel('Z [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
title('Tip Orientation Over Time', 'Fontname', 'Times New Roman', 'fontsize', 16);
hold on;

interval = 50;

% Original data (Yu)
posX_orig = Yu{1,2}(2,1:interval:end);
posY_orig = Yu{1,2}(3,1:interval:end);
posZ_orig = Yu{1,2}(1,1:interval:end);
quate1_orig = Yu{1,2}(4,1:interval:end);
quate2_orig = Yu{1,2}(5,1:interval:end);
quate3_orig = Yu{1,2}(6,1:interval:end);
quate4_orig = Yu{1,2}(7,1:interval:end);

% Remapped data (yDataRemapped)
posX_est = yDataRemapped{1,2}(2,1:interval:end);
posY_est = yDataRemapped{1,2}(3,1:interval:end);
posZ_est = yDataRemapped{1,2}(1,1:interval:end);
quate1_est = yDataRemapped{1,2}(4,1:interval:end);
quate2_est = yDataRemapped{1,2}(5,1:interval:end);
quate3_est = yDataRemapped{1,2}(6,1:interval:end);
quate4_est = yDataRemapped{1,2}(7,1:interval:end);

% Calculate x-axis direction vector for original data
x_axis_world_orig = zeros(length(quate1_orig), 3);
for i = 1:length(quate1_orig)
    q = [quate1_orig(i), quate2_orig(i), quate3_orig(i), quate4_orig(i)]; % [w x y z]
    % Quaternion rotation matrix
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),   2*(x*z+w*y);
         2*(x*y+w*z), 1-2*(x^2+z^2),   2*(y*z-w*x);
         2*(x*z-w*y),   2*(y*z+w*x), 1-2*(x^2+y^2)];
    % X-axis: Points along the length of the sensor body — usually from the cable end toward the tip.
    x_axis_world_orig(i,:) = R * [1;0;0];
end

% Calculate x-axis direction vector for remapped data
x_axis_world_est = zeros(length(quate1_est), 3);
for i = 1:length(quate1_est)
    q = [quate1_est(i), quate2_est(i), quate3_est(i), quate4_est(i)]; % [w x y z]
    % Quaternion rotation matrix
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),   2*(x*z+w*y);
         2*(x*y+w*z), 1-2*(x^2+z^2),   2*(y*z-w*x);
         2*(x*z-w*y),   2*(y*z+w*x), 1-2*(x^2+y^2)];
    % X-axis: Points along the length of the sensor body — usually from the cable end toward the tip.
    x_axis_world_est(i,:) = R * [1;0;0];
end

scale = 1;

% Plot original data arrows (Yu) - Red/Orange color
h1_arrow = quiver3(posX_orig(1), posY_orig(1), -posZ_orig(1), ...
        x_axis_world_orig(1,2), x_axis_world_orig(1,3), -x_axis_world_orig(1,1), ...
        scale, 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'MaxHeadSize', 0.8);
for i = 2:length(posX_orig)
    quiver3(posX_orig(i), posY_orig(i), -posZ_orig(i), ...
            x_axis_world_orig(i,2), x_axis_world_orig(i,3), -x_axis_world_orig(i,1), ...
            scale, 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'MaxHeadSize', 0.8);
end

% Plot remapped data arrows (yDataRemapped) - Blue color with offset
offset = 0.1; % Small offset to avoid complete overlap
h2_arrow = quiver3(posX_est(1)+offset, posY_est(1)+offset, -posZ_est(1), ...
        x_axis_world_est(1,2), x_axis_world_est(1,3), -x_axis_world_est(1,1), ...
        scale*0.8, 'Color', [0.1 0.3 0.8], 'LineWidth', 1.5, 'MaxHeadSize', 1, 'LineStyle', '--');
for i = 2:length(posX_est)
    quiver3(posX_est(i)+offset, posY_est(i)+offset, -posZ_est(i), ...
            x_axis_world_est(i,2), x_axis_world_est(i,3), -x_axis_world_est(i,1), ...
            scale*0.8, 'Color', [0.1 0.3 0.8], 'LineWidth', 1.5, 'MaxHeadSize', 1, 'LineStyle', '--');
end

% Add legend
legend([h1_arrow, h2_arrow], {'Original Poses', 'Reconstructed Poses'}, ...
       'Location', 'best', 'FontSize', 14);

zlim([-1,1]);
view([-30 40]);
set(gca,'linewidth',2);
ax = gca;
ax.GridLineWidth = 2;
set(gca,'FontSize',16);
grid on;
axis equal;