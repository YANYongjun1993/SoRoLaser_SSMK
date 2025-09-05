clearvars;
clc
clear all
close all
%% train and validation dataset frequency alalysis
folder = 'E:\Research\SRIL\14th\TroRevision\Figure\OpenLoopComp\smallData\SSMpoly\processedData';
% Load lowDimCtrl
fullFileName = fullfile(folder, 'softKpmRedu_train-14_val-3_2025-06-18.mat');
tmp = load(fullFileName);
if isfield(tmp, 'val')
    validationDataSmall = tmp.val;
else
    error('Variable val not found in %s', fullFileName);
end
if isfield(tmp, 'train')
    trainDataSmall = tmp.train;
else
    error('Variable train not found in %s', fullFileName);
end

folder = 'E:\Research\SRIL\14th\TroRevision\Figure\OpenLoopComp\bigData\Koopman\processedData';
% Load lowDimCtrl
fullFileName = fullfile(folder, 'softKpmObs_capacity_train-17_val-3_2025-06-25.mat');
tmp = load(fullFileName);
if isfield(tmp, 'val')
    validationDataBig = tmp.val;
else
    error('Variable val not found in %s', fullFileName);
end
if isfield(tmp, 'train')
    trainDataBig = tmp.train;
else
    error('Variable train not found in %s', fullFileName);
end
%% Plot FFT results of trainDataBig
fs = 100; % Sampling frequency in Hz (change if known)
nTrain = size(trainDataBig, 2);

figure;
width=350; % width of the figure, unit is pixel
height=250; % height of the figure, unit is pixel
left=50;%
bottem=50;%
set(gcf,'position',[left,bottem,width,height]);
hold on;
color1 = [0 0.4470 0.7410]; % blue for 1:13
color2 = [0.8500 0.3250 0.0980]; % orange for 14:17

for i = 1:nTrain
    if i <= 13
        thisColor = color1;
    else
        thisColor = color2;
    end
    y = trainDataBig{1, i}.y(:,1); % Only plot y(:,1)
    N = length(y);
    f = (0:N-1)*(fs/N); % Frequency vector

    Y = fft(y);
    P2 = abs(Y/N); % Two-sided spectrum
    P1 = P2(1:floor(N/2)+1); % Single-sided
    P1(2:end-1) = 2*P1(2:end-1);
    f_plot = f(1:floor(N/2)+1);

    % Only plot up to 5 Hz
    idx = f_plot <= 5;
    plot(f_plot(idx), P1(idx), 'Color', thisColor, 'LineWidth', 1);
end
hold off;
xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',12);
ylabel('Amplitude','Fontname','Times New Roman','fontsize',12);
title('FFT of New Training Dataset','Fontname','Times New Roman','fontsize',12);
set(gca,'linewidth',1);
ax = gca;
ax.GridLineWidth = 1;
set(gca,'FontSize',12);
xlim([0 5]);
ylim([0 0.15]);
grid on;
box on;
% legend({'Autonomous Dynamics','Controlled Trajectory'}, 'Location', 'northeast','Fontname','Times New Roman','fontsize',12);
%% Plot FFT results of validationDataBig
fs = 100; % Sampling frequency in Hz (change if known)
nVal = size(validationDataBig, 2);

figure;
width=350; % width of the figure, unit is pixel
height=250; % height of the figure, unit is pixel
left=50;%
bottem=50;%
set(gcf,'position',[left,bottem,width,height]);
hold on;
color1 = [0 0.4470 0.7410]; % blue for 1:13
color2 = [0.8500 0.3250 0.0980]; % orange for 14:17

for i = 1:nVal
    if i <= 2
        thisColor = color1;
    else
        thisColor = color2;
    end
    y = validationDataBig{1, i}.y(:,1); % Only plot y(:,1)
    N = length(y);
    f = (0:N-1)*(fs/N); % Frequency vector

    Y = fft(y);
    P2 = abs(Y/N); % Two-sided spectrum
    P1 = P2(1:floor(N/2)+1); % Single-sided
    P1(2:end-1) = 2*P1(2:end-1);
    f_plot = f(1:floor(N/2)+1);

    % Only plot up to 5 Hz
    idx = f_plot <= 5;
    plot(f_plot(idx), P1(idx), 'Color', thisColor, 'LineWidth', 1);
end
hold off;
xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',12);
ylabel('Amplitude','Fontname','Times New Roman','fontsize',12);
title('FFT of New Validation Dataset','Fontname','Times New Roman','fontsize',12);
set(gca,'linewidth',1);
ax = gca;
ax.GridLineWidth = 1;
set(gca,'FontSize',12);
xlim([0 5]);
ylim([0 0.15]);
grid on;
box on;
%% Plot FFT results of trainDataSmall
fs = 100; % Sampling frequency in Hz (change if known)
nTrain = size(trainDataSmall, 2);

figure;
width=350; % width of the figure, unit is pixel
height=250; % height of the figure, unit is pixel
left=50;%
bottem=50;%
set(gcf,'position',[left,bottem,width,height]);
hold on;
color1 = [0 0.4470 0.7410]; % blue for 1:13
color2 = [0.8500 0.3250 0.0980]; % orange for 14:17

for i = 1:nTrain
    if i <= 13
        thisColor = color1;
    else
        thisColor = color2;
    end
    y = trainDataSmall{1, i}.x(:,1); % Only plot y(:,1)
    N = length(y);
    f = (0:N-1)*(fs/N); % Frequency vector

    Y = fft(y);
    P2 = abs(Y/N); % Two-sided spectrum
    P1 = P2(1:floor(N/2)+1); % Single-sided
    P1(2:end-1) = 2*P1(2:end-1);
    f_plot = f(1:floor(N/2)+1);

    % Only plot up to 5 Hz
    idx = f_plot <= 5;
    plot(f_plot(idx), P1(idx), 'Color', thisColor, 'LineWidth', 1);
end
hold off;
xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',12);
ylabel('Amplitude','Fontname','Times New Roman','fontsize',12);
% title('FFT of Training Dataset','Fontname','Times New Roman','fontsize',12);
set(gca,'linewidth',1);
ax = gca;
ax.GridLineWidth = 1;
set(gca,'FontSize',12);
xlim([0 5]);
ylim([0 0.15]);
grid on;
box on;
%% Plot FFT results of validationDataSmall
fs = 100; % Sampling frequency in Hz (change if known)
nVal = size(validationDataSmall, 2);

figure;
width=350; % width of the figure, unit is pixel
height=250; % height of the figure, unit is pixel
left=50;%
bottem=50;%
set(gcf,'position',[left,bottem,width,height]);
hold on;
color1 = [0 0.4470 0.7410]; % blue for 1:13
color2 = [0.8500 0.3250 0.0980]; % orange for 14:17

for i = 1:nVal
    if i <= 2
        thisColor = color1;
    else
        thisColor = color2;
    end
    y = validationDataSmall{1, i}.x(:,1); % Only plot y(:,1)
    N = length(y);
    f = (0:N-1)*(fs/N); % Frequency vector

    Y = fft(y);
    P2 = abs(Y/N); % Two-sided spectrum
    P1 = P2(1:floor(N/2)+1); % Single-sided
    P1(2:end-1) = 2*P1(2:end-1);
    f_plot = f(1:floor(N/2)+1);

    % Only plot up to 5 Hz
    idx = f_plot <= 5;
    plot(f_plot(idx), P1(idx), 'Color', thisColor, 'LineWidth', 1);
end
hold off;
xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',12);
ylabel('Amplitude','Fontname','Times New Roman','fontsize',12);
% title('FFT of New Validation Dataset','Fontname','Times New Roman','fontsize',12);
set(gca,'linewidth',1);
ax = gca;
ax.GridLineWidth = 1;
set(gca,'FontSize',12);
xlim([0 5]);
ylim([0 0.15]);
grid on;
box on;