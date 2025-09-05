%% Create figure for orientation
clear all
clc
load Figure8BasicSpdRegFactor1.mat

%% Delay-embedding
SSMDim = 2;
overEmbed = 32;
[figure8DataDelay, opts_embd] = coordinatesEmbeddingControl(zu, SSMDim, 'OverEmbedding', overEmbed);
embedDim = size(figure8DataDelay{1,2},1);
outdofsDelay = [embedDim-6 embedDim-5 embedDim-4 embedDim-3 embedDim-2 embedDim-1 embedDim];
%%
load SSM_model.mat
CtrlData1 = projectTrajectories(ssm_model.IMInfoCtrl, zu); % Project down
obsRec1 = ssm_model.IMInfoCtrl.parametrization.map(CtrlData1{1,2});
lowDimDelay = projectTrajectories(ssm_model.IMInfo, figure8DataDelay); % Project down
newLowDim = ssm_model.IMInfo.trans.map(lowDimDelay{1,2});
obsRec2 = ssm_model.IMInfoCtrl.parametrization.map(newLowDim);
slice_duration = size(obsRec2,2);
shift = 10;
slice_begain = size(zu{1,1},2) - slice_duration + 1 - shift;
slice_end = size(zu{1,2},2)-shift;
%% Extract positions and quaternions
posX_est_1 = zu{1,2}(2,slice_begain:slice_end);
posY_est_1 = zu{1,2}(3,slice_begain:slice_end);
posZ_est_1 = zu{1,2}(1,slice_begain:slice_end);
time_steps_1 = zu{1,1}(1,1:slice_duration) - zu{1,1}(1,1);
posX_remap_1 = obsRec2(2,:);
posY_remap_1 = obsRec2(3,:);
posZ_remap_1 = obsRec2(1,:);
%%
load Figure8BasicSpdRegFactor10.mat

%% Delay-embedding
SSMDim = 2;
overEmbed = 32;
[figure8DataDelay, opts_embd] = coordinatesEmbeddingControl(zu, SSMDim, 'OverEmbedding', overEmbed);
embedDim = size(figure8DataDelay{1,2},1);
outdofsDelay = [embedDim-6 embedDim-5 embedDim-4 embedDim-3 embedDim-2 embedDim-1 embedDim];

%%
load SSM_model.mat
CtrlData1 = projectTrajectories(ssm_model.IMInfoCtrl, zu); % Project down
obsRec1 = ssm_model.IMInfoCtrl.parametrization.map(CtrlData1{1,2});
lowDimDelay = projectTrajectories(ssm_model.IMInfo, figure8DataDelay); % Project down
newLowDim = ssm_model.IMInfo.trans.map(lowDimDelay{1,2});
obsRec2 = ssm_model.IMInfoCtrl.parametrization.map(newLowDim);
slice_duration = size(obsRec2,2);
shift = 10;
slice_begain = size(zu{1,1},2) - slice_duration + 1 - shift;
slice_end = size(zu{1,2},2)-shift;
%% Extract positions and quaternions
posX_est_10 = zu{1,2}(2,slice_begain:slice_end);
posY_est_10 = zu{1,2}(3,slice_begain:slice_end);
posZ_est_10 = zu{1,2}(1,slice_begain:slice_end);
time_steps_10 = zu{1,1}(1,1:slice_duration) - zu{1,1}(1,1);
posX_remap_10 = obsRec2(2,:);
posY_remap_10 = obsRec2(3,:);
posZ_remap_10 = obsRec2(1,:);
%%
load Figure8BasicSpdRegFactor15.mat
%% Delay-embedding
SSMDim = 2;
overEmbed = 32;
[figure8DataDelay, opts_embd] = coordinatesEmbeddingControl(zu, SSMDim, 'OverEmbedding', overEmbed);
embedDim = size(figure8DataDelay{1,2},1);
outdofsDelay = [embedDim-6 embedDim-5 embedDim-4 embedDim-3 embedDim-2 embedDim-1 embedDim];
%%
load SSM_model.mat
CtrlData1 = projectTrajectories(ssm_model.IMInfoCtrl, zu); % Project down
obsRec1 = ssm_model.IMInfoCtrl.parametrization.map(CtrlData1{1,2});
lowDimDelay = projectTrajectories(ssm_model.IMInfo, figure8DataDelay); % Project down
newLowDim = ssm_model.IMInfo.trans.map(lowDimDelay{1,2});
obsRec2 = ssm_model.IMInfoCtrl.parametrization.map(newLowDim);
slice_duration = size(obsRec2,2);
shift = 10;
slice_begain = size(zu{1,1},2) - slice_duration + 1 - shift;
slice_end = size(zu{1,2},2)-shift;
%% Extract positions and quaternions
posX_est_15 = zu{1,2}(2,slice_begain:slice_end);
posY_est_15 = zu{1,2}(3,slice_begain:slice_end);
posZ_est_15 = zu{1,2}(1,slice_begain:slice_end);
time_steps_15 = zu{1,1}(1,1:slice_duration) - zu{1,1}(1,1);
posX_remap_15 = obsRec2(2,:);
posY_remap_15 = obsRec2(3,:);
posZ_remap_15 = obsRec2(1,:);
%%
figure('Position', [100, 100, 700, 700]);
subplot(2,1,1)
plot(time_steps_1,posX_est_1,'-','color',[0.4 0.6 0.9],'LineWidth',2)
hold on
plot(time_steps_1,posX_remap_1,'--','color',[0.1 0.3 0.8],'LineWidth',3)
hold on
plot(time_steps_10,posX_est_10,'-','color',[1 0.7 0.4],'LineWidth',2)
hold on
plot(time_steps_10,posX_remap_10,'--','color',[1 0.5 0],'LineWidth',2)
hold on
plot(time_steps_15,posX_est_15,'-','color',[1 0.4 0.4],'LineWidth',2)
hold on
plot(time_steps_15,posX_remap_15,'--','color',[0.8 0.1 0.1],'LineWidth',2)
% xlabel('Time [s]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
ylabel('X [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
set(gca,'linewidth',2);
ax = gca;
ax.GridLineWidth = 2;
set(gca,'FontSize',16);
grid on;
subplot(2,1,2)
plot(time_steps_1,posY_est_1,'-','color',[0.4 0.6 0.9],'LineWidth',2)
hold on
plot(time_steps_1,posY_remap_1,'--','color',[0.1 0.3 0.8],'LineWidth',2)
hold on
plot(time_steps_10,posY_est_10,'-','color',[1 0.7 0.4],'LineWidth',2)
hold on
plot(time_steps_10,posY_remap_10,'--','color',[1 0.5 0],'LineWidth',2)
hold on
plot(time_steps_15,posY_est_15,'-','color',[1 0.4 0.4],'LineWidth',2)
hold on
plot(time_steps_15,posY_remap_15,'--','color',[0.8 0.1 0.1],'LineWidth',2)
xlabel('Time [s]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Y [mm]', 'Fontname', 'Times New Roman', 'fontsize', 16, 'fontweight', 'bold');
legend("Durations:5 s — EM Sensor (Original)","Durations:5 s — SSM (Reconstructed)","Durations:2.5 s — EM Sensor (Original)","Durations:2.5 s — SSM (Reconstructed)","Durations:2 s — EM Sensor (Original)","Durations:2 s — SSM (Reconstructed)", 'FontName', 'Times New Roman', 'NumColumns', 2)
set(gca,'linewidth',2);
ax = gca;
ax.GridLineWidth = 2;
set(gca,'FontSize',16);
grid on;
%% Concatenating the trajectories
pos_est_all = [posX_est_1,posX_est_10,posX_est_15,posY_est_1,posY_est_10,posY_est_15,posZ_est_1,posZ_est_10,posZ_est_15];
pos_remap_all = [posX_remap_1,posX_remap_10,posX_remap_15,posY_remap_1,posY_remap_10,posY_remap_15,posZ_remap_1,posZ_remap_10,posZ_remap_15];
% pos_est_all = [posX_est_1,posX_est_10,posX_est_15,posY_est_1,posY_est_10,posY_est_15];
% pos_remap_all = [posX_remap_1,posX_remap_10,posX_remap_15,posY_remap_1,posY_remap_10,posY_remap_15];
load dominant.mat
% decayComb = [];
% for i = 1:size(zDataDelay,1)
%     decayComb = [decayComb,zDataDelay{i,2}];
% end
decayCombFinal = [zDataDelay{1,2}(1,:),zDataDelay{1,2}(2,:)];
%% Plot FFT results of trainDataBig
fs = 100; % Sampling frequency in Hz (change if known)
fs_delay = 1000;
figure('Position', [100, 100, 700, 400]);
hold on;
color1 = [0 0.4470 0.7410]; % blue
color2 = [0.8500 0.3250 0.0980]; % orange
%%
y = pos_est_all; % Only plot y(:,1)
N = length(y);
f = (0:N-1)*(fs/N); % Frequency vector
Y = fft(y);
P2 = abs(Y/N); % Two-sided spectrum
P1 = P2(1:floor(N/2)+1); % Single-sided
P1(2:end-1) = 2*P1(2:end-1);
f_plot = f(1:floor(N/2)+1);

% Only plot up to 5 Hz
idx = f_plot <= 50;
plot(f_plot(idx), P1(idx), 'Color', color1, 'LineWidth', 1);
%%
y = pos_remap_all; % Only plot y(:,1)
N = length(y);
f = (0:N-1)*(fs/N); % Frequency vector
Y = fft(y);
P2 = abs(Y/N); % Two-sided spectrum
P1 = P2(1:floor(N/2)+1); % Single-sided
P1(2:end-1) = 2*P1(2:end-1);
f_plot = f(1:floor(N/2)+1);

% Only plot up to 5 Hz
idx = f_plot <= 50;
plot(f_plot(idx), P1(idx), 'Color', color2, 'LineWidth', 1);
%%
xline(3.1, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

hold off;
xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',16, 'fontweight', 'bold');
ylabel('Amplitude','Fontname','Times New Roman','fontsize',16, 'fontweight', 'bold');
title('FFT of New Training Dataset','Fontname','Times New Roman','fontsize',12, 'fontweight', 'bold');
legend("EM Sensor (Original)"," SSM (Reconstructed)", 'FontName', 'Times New Roman')

set(gca,'linewidth',2);
ax = gca;
ax.GridLineWidth = 2;
set(gca,'FontSize',16);
% xlim([0 5]);
ylim([0 0.15]);
grid on;
box on;
% legend({'Autonomous Dynamics','Controlled Trajectory'}, 'Location', 'northeast','Fontname','Times New Roman','fontsize',12);