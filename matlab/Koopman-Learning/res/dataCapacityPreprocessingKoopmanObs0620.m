clearvars;
clc
clear all
close all
%% decay data--autunomous dynamics
load rawData/decaydata202503271424.mat % with hand
sensor1_Quater1 = get(logsout,"sensor1_Quater1");
time = sensor1_Quater1.Values.Time;
Sensor1Quater1 = sensor1_Quater1.Values.data;
sensor1_Quater2 = get(logsout,"sensor1_Quater2");
Sensor1Quater2 = sensor1_Quater2.Values.data;
sensor1_Quater3 = get(logsout,"sensor1_Quater3");
Sensor1Quater3 = sensor1_Quater3.Values.data;
sensor1_Quater4 = get(logsout,"sensor1_Quater4");
Sensor1Quater4 = sensor1_Quater4.Values.data;
sensor1_x = get(logsout,"sensor1_x");
Sensor1X = sensor1_x.Values.data;
sensor1_y = get(logsout,"sensor1_y");
Sensor1Y = sensor1_y.Values.data;
sensor1_z = get(logsout,"sensor1_z");
Sensor1Z = sensor1_z.Values.data;
sensor2_Quater1 = get(logsout,"sensor2_Quater1");
Sensor2Quater1 = sensor2_Quater1.Values.data;
sensor2_Quater2 = get(logsout,"sensor2_Quater2");
Sensor2Quater2 = sensor2_Quater2.Values.data;
sensor2_Quater3 = get(logsout,"sensor2_Quater3");
Sensor2Quater3 = sensor2_Quater3.Values.data;
sensor2_Quater4 = get(logsout,"sensor2_Quater4");
Sensor2Quater4 = sensor2_Quater4.Values.data;
sensor2_x = get(logsout,"sensor2_x");
Sensor2X = sensor2_x.Values.data;
sensor2_y = get(logsout,"sensor2_y");
Sensor2Y = sensor2_y.Values.data;
sensor2_z = get(logsout,"sensor2_z");
Sensor2Z = sensor2_z.Values.data;
% concatenating the position and quaterion
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];
%% Slice the whole data as some segments, "z" represents reduced coordinate, "o" represent observable coordinate.
% 2025:0604:1439 reextract the decaying data to make the range of the position under 20 mm
nInitPoints = 15;
samplingInterval = 0.001;
oData= cell(nInitPoints,2);
slice1 = [8985,10341];
oData{1,1} = 0:samplingInterval:(slice1(2) - slice1(1))*samplingInterval;
% oData11 = 0:samplingInterval:(slice1(2) - slice1(1))*samplingInterval;
oData{1,2} = [Sensor1(:,slice1(1):slice1(2))];
% oData12 = [pos_SoMlp_1(1,slice1(1):slice1(2))];
slice2 = [13189,14318];
oData{2,1} = 0:samplingInterval:(slice2(2) - slice2(1))*samplingInterval;
oData{2,2} = [Sensor1(:,slice2(1):slice2(2))];
slice3 = [17251, 18677];
oData{3,1} = 0:samplingInterval:(slice3(2) - slice3(1))*samplingInterval;
oData{3,2} = [Sensor1(:,slice3(1):slice3(2))];
slice4 = [21682, 23300];
oData{4,1} = 0:samplingInterval:(slice4(2) - slice4(1))*samplingInterval;
oData{4,2} = [Sensor1(:,slice4(1):slice4(2))];
slice5 = [26144, 27281];
oData{5,1} = 0:samplingInterval:(slice5(2) - slice5(1))*samplingInterval;
oData{5,2} = [Sensor1(:,slice5(1):slice5(2))];
slice6 = [30621, 31690];
oData{6,1} = 0:samplingInterval:(slice6(2) - slice6(1))*samplingInterval;
oData{6,2} = [Sensor1(:,slice6(1):slice6(2))];
slice7 = [34722, 36318];
oData{7,1} = 0:samplingInterval:(slice7(2) - slice7(1))*samplingInterval;
oData{7,2} = [Sensor1(:,slice7(1):slice7(2))];
slice8 = [38717, 40577];
oData{8,1} = 0:samplingInterval:(slice8(2) - slice8(1))*samplingInterval;
oData{8,2} = [Sensor1(:,slice8(1):slice8(2))];
slice9 = [42841, 44114];
oData{9,1} = 0:samplingInterval:(slice9(2) - slice9(1))*samplingInterval;
oData{9,2} = [Sensor1(:,slice9(1):slice9(2))];
slice10 = [47185, 48151];
oData{10,1} = 0:samplingInterval:(slice10(2) - slice10(1))*samplingInterval;
oData{10,2} = [Sensor1(:,slice10(1):slice10(2))];
slice11 = [50683, 51695];
oData{11,1} = 0:samplingInterval:(slice11(2) - slice11(1))*samplingInterval;
oData{11,2} = [Sensor1(:,slice11(1):slice11(2))];
slice12 = [54580, 55975];
oData{12,1} = 0:samplingInterval:(slice12(2) - slice12(1))*samplingInterval;
oData{12,2} = [Sensor1(:,slice12(1):slice12(2))];
slice13 = [58369, 59398];
oData{13,1} = 0:samplingInterval:(slice13(2) - slice13(1))*samplingInterval;
oData{13,2} = [Sensor1(:,slice13(1):slice13(2))];
slice14 = [61526, 62897];
oData{14,1} = 0:samplingInterval:(slice12(2) - slice12(1))*samplingInterval;
oData{14,2} = [Sensor1(:,slice12(1):slice12(2))];
slice15 = [64795, 65569];
oData{15,1} = 0:samplingInterval:(slice13(2) - slice13(1))*samplingInterval;
oData{15,2} = [Sensor1(:,slice13(1):slice13(2))];
%% Truncate and transform trajectories to the origin
nPoints = size(oData,1);
oDataTrunc= cell(nPoints,2);

for i = 1:nPoints
    oDataTrunc{i,1} = oData{i,1};
    oDataTrunc{i,2} = oData{i,2};
    for iFeature = 1:size(oDataTrunc{i,2},1)
        oDataTrunc{i,2}(iFeature,:) = oDataTrunc{i,2}(iFeature,:) - oDataTrunc{i,2}(iFeature,end);
    end
end
%% Multi-stage smoothing for reliable differential calculus
% Parameters
windowSize = 51;  % Must be odd
polynomialOrder = 2;  % Balance between smoothing and feature preservation
medianFilterSize = 5;  % For removing spikes
butterworthOrder = 4;  % Butterworth filter order
cutoffFreq = 0.1;     % Normalized cutoff frequency (0 to 1)

oDataSmooth = cell(nPoints, 2);

for iTraj = 1:nPoints
    % Copy time data
    oDataSmooth{iTraj,1} = oDataTrunc{iTraj,1};
    
    % Get position data
    pos_data = oDataTrunc{iTraj,2};
    smoothed_pos = zeros(size(pos_data));
    
    for iDim = 1:size(pos_data,1)
        % 1. Apply median filter to remove spikes
        temp_data = medfilt1(pos_data(iDim,:), medianFilterSize);
        
        % 2. Apply Butterworth filter
        [b, a] = butter(butterworthOrder, cutoffFreq, 'low');
        temp_data = filtfilt(b, a, temp_data);
        
        % 3. Final smoothing with Savitzky-Golay
        smoothed_pos(iDim,:) = sgolayfilt(temp_data, polynomialOrder, windowSize);
    end
    
    % Store smoothed data
    oDataSmooth{iTraj,2} = smoothed_pos;
end
%% Delay-embedding
SSMDim = 4;
overEmbed = 45;
[yData, opts_embd] = coordinatesEmbeddingControl(oDataSmooth, SSMDim, 'OverEmbedding', overEmbed);
embedDim = size(yData{1,2},1);
outdofsDelay = [embedDim-6 embedDim-5 embedDim-4 embedDim-3 embedDim-2 embedDim-1 embedDim];
oDataDelay = yData;
%% Show PCA (principal component analysis) infos
disp('PCA')
oSnapshots = cat(2,oDataDelay{:,2});
% Perform SVD on displacement field
[V,S,U] = svd(oSnapshots, 0);
% Note we assume data centered around the origin, which is the
% fixed point of our system
customFigure;
maxModesPlot = 10;
l2vals = diag(S).^2;
plot(1:length(l2vals),cumsum(l2vals)/sum(l2vals)*100,'.-','Linewidth',2,'MarkerSize',12)
xlabel('Number of Modes');
ylabel('Variance [%]');
xlim([1 maxModesPlot]);
set(gca,'XTick',1:1:maxModesPlot);
Vde = V(:,1:SSMDim);
% Set reduced-coordinates as PCA projection to the first four modes
zDataDelay = oDataDelay;
for iTraj = 1:nPoints
    zDataDelay{iTraj,2} = transpose(Vde)*oDataDelay{iTraj,2};
end 
%% Set test and train trajectories out of all the trajectories
indOut = [3 9 13 15];
indTest = [2 8 11];
indTrain = setdiff(1:size(oData, 1),union(indTest,indOut));
%% Geometry: Manifold fitting with prescribed graph. Higher orders tend to overfit.
SSMOrder = 1;
obsDim = 6;
Nfolds = 5;
Nregvals = 30; 
if SSMOrder == 1
    % Need to recompute as the columns of Vde might change in order
    IMInfo = IMGeometry(oDataDelay(indTrain,:), SSMDim, SSMOrder);
    Vde = IMInfo.parametrization.tangentSpaceAtOrigin;
    zDataDelay = oDataDelay;
    for iTraj = 1:size(oData, 1)
        zDataDelay{iTraj,2} = transpose(Vde)*oDataDelay{iTraj,2};
    end
else
    % IMInfo = IMGeometryCV(oDataTrunc, SSMDim, SSMOrder,'reducedCoordinates',zDataDelay,'R_PolyOrd', ROMOrder, 'n_folds', Nfolds, 'l_vals', logspace(-6,0,Nregvals), 'style', 'normalform');
    % IMInfoInv = IMGeometry(etaDataTrunc,obsDim, SSMOrder,'reducedCoordinates',zDataDelay);
    IMInfo = IMGeometry(oDataDelay(indTrain,:), SSMDim, SSMOrder,'reducedCoordinates',zDataDelay(indTrain,:));
    % IMInfoInv = IMGeometry(zDataDelay,obsDim, SSMOrder,'reducedCoordinates',tipDataCrt);
    IMInfo.chart.map = @(x) transpose(Vde)*x;
end
%% Switch to observable coordinates
oDataCrtObs = oDataDelay;
for iTraj = 1:size(yData,1)
    oDataCrtObs{iTraj,2} = oDataDelay{iTraj,2}(outdofsDelay,:);
end
IMInfo_paramonly = IMGeometry(oDataCrtObs(indTrain,:), SSMDim, SSMOrder,'reducedCoordinates',zDataDelay(indTrain,:),'style','custom');
tanSpace0_not_orth = IMInfo_paramonly.parametrization.tangentSpaceAtOrigin;
tanSpace0 = orthogonalizeGramSchmidt(tanSpace0_not_orth);
% Change coordinates
zDataCrt = zDataDelay;
for iTraj = 1:size(zDataDelay,1)
    zDataCrt{iTraj,2} = transpose(tanSpace0)*tanSpace0_not_orth*zDataDelay{iTraj,2};
end
% Fit new parametrization
IMInfo_paramonly = IMGeometry(oDataCrtObs(indTrain,:), SSMDim, SSMOrder,'reducedCoordinates',zDataCrt(indTrain,:),'Ve',tanSpace0,'style','custom');
disp('Results of new parametrization fit:')
oRec = liftTrajectories(IMInfo_paramonly, zDataCrt);
normedTrajDist = computeTrajectoryErrors(oRec, oDataCrtObs);
meanErrorGeoTrain = mean(normedTrajDist(indTrain))*100
meanErrorGeoTest = mean(normedTrajDist(indTest))*100

IMInfoCtrl.chart.map = @(x) transpose(tanSpace0)*x;
IMInfoCtrl = struct('chart', IMInfoCtrl.chart, 'parametrization', IMInfo_paramonly.parametrization);
% repack the data
dataDecay = cell(1,nPoints);
for i = 1:nPoints
    % Assign struct to the i-th cell
    dataDecay{i} = struct('t', zDataCrt{i,1}', 'x', zDataCrt{i,2}', 'y', oDataDelay{i,2}(outdofsDelay,:)');
end
%% Export mapping coefficients and input matrix
ssm_model = struct();
ssm_model.sizeX = 4;
ssm_model.sizeU = 4;
ssm_model.SSMDim = SSMDim;
ssm_model.IMInfoCtrl = IMInfoCtrl;
save('SSM_model.mat', 'ssm_model', '-v7');
% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('processedData', sprintf('koopmanDecayObs_%s.mat', current_date));

% Save the trajectory data
save(fileName, 'dataDecay');
%% load new control data 1
clearvars -except dataDecay;
clc
close all
% Load new controller data with random inputs
load rawData/202505252158.mat

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
load SSM_model.mat
CtrlData1 = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
% pack the trajectory and control input data at the low dimentional as struct
CtrlData1Packed = struct('t', zuSmooth{1,1}', 'u', uCell{1,2}', 'x', CtrlData1{1,2}','y',Yu{1,2}');
% save('Yu202505252158.mat', 'CtrlData1Packed');
%% load new control data 2
clearvars -except dataDecay CtrlData1Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202505201626.mat
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
%% setup matrix
zu = {};
keptIdxs = 5429:25555;
y_eq = Sensor2(:,5000);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor2(:, keptIdxs) - y_eq);

uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorPos(:, keptIdxs); 
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
load SSM_model.mat
CtrlData2 = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
% pack the trajectory and control input data at the low dimentional as struct
CtrlData2Packed = struct('t', zuSmooth{1,1}', 'u', uCell{1,2}',  'x', CtrlData2{1,2}','y',Yu{1,2}');
%% load new control data 3
clearvars -except dataDecay CtrlData1Packed CtrlData2Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202505252202.mat
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
%% setup matrix
zu = {};
keptIdxs = 5487:15392;
y_eq = Sensor2(:,15600);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor2(:, keptIdxs) - y_eq);

uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorPos(:, keptIdxs); 
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
load SSM_model.mat
CtrlData3 = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
% pack the trajectory and control input data at the low dimentional as struct
CtrlData3Packed = struct('t', zuSmooth{1,1}', 'u', uCell{1,2}',  'x', CtrlData3{1,2}', 'y',Yu{1,2}');
%% load new control data 1:trajectory_data_20250611.mat, Gain=800, spdRegFactor=1
clearvars -except dataDecay CtrlData1Packed CtrlData2Packed CtrlData3Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202506112038.mat

sensor1_Quater1 = get(logsout,"sensor1_Quater1");
time = sensor1_Quater1.Values.Time;
Sensor1Quater1 = sensor1_Quater1.Values.data;
sensor1_Quater2 = get(logsout,"sensor1_Quater2");
Sensor1Quater2 = sensor1_Quater2.Values.data;
sensor1_Quater3 = get(logsout,"sensor1_Quater3");
Sensor1Quater3 = sensor1_Quater3.Values.data;
sensor1_Quater4 = get(logsout,"sensor1_Quater4");
Sensor1Quater4 = sensor1_Quater4.Values.data;
sensor1_x = get(logsout,"sensor1_x");
Sensor1X = sensor1_x.Values.data;
sensor1_y = get(logsout,"sensor1_y");
Sensor1Y = sensor1_y.Values.data;
sensor1_z = get(logsout,"sensor1_z");
Sensor1Z = sensor1_z.Values.data;
% concatenating the position and quaterion
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];

flag = get(logsout,"flag");
flagStatus = flag.Values.data;

Motor1_pos = get(logsout,"Motor1_pos_req_qc_mode2");
Motor1Pos = Motor1_pos.Values.data;
Motor2_pos = get(logsout,"Motor2_pos_req_qc_mode2");
Motor2Pos = Motor2_pos.Values.data;
Motor3_pos = get(logsout,"Motor3_pos_req_qc_mode2");
Motor3Pos = Motor3_pos.Values.data;
Motor4_pos = get(logsout,"Motor4_pos_req_qc_mode2");
Motor4Pos = Motor4_pos.Values.data;
uMotorPos = [Motor1Pos';Motor2Pos';Motor3Pos';Motor4Pos'];
%% setup matrix
zu = {};
keptIdxs = 5570:15562;
y_eq = Sensor1(:,15800);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);
u_eq = uMotorPos(:, 15800);
uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorPos(:, keptIdxs) - u_eq; 
%% Multi-stage smoothing for reliable differential calculus
% Parameters
windowSize = 101;  % Must be odd
polynomialOrder = 2;  % Balance between smoothing and feature preservation
medianFilterSize = 7;  % For removing spikes
butterworthOrder = 4;  % Butterworth filter order
cutoffFreq = 0.01;     % Normalized cutoff frequency (0 to 1)
  
for iDim = 1:size(Sensor1,1)
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
load SSM_model.mat
CtrlData4 = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
% pack the trajectory and control input data at the low dimentional as struct
CtrlData4Packed = struct('t', zuSmooth{1,1}', 'u', uCell{1,2}', 'x', CtrlData4{1,2}','y',Yu{1,2}');
%% load new control data 2:trajectory_data_20250611.mat, Gain=800, spdRegFactor=1.5
clearvars -except dataDecay CtrlData1Packed CtrlData2Packed CtrlData3Packed CtrlData4Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202506112044.mat

sensor1_Quater1 = get(logsout,"sensor1_Quater1");
time = sensor1_Quater1.Values.Time;
Sensor1Quater1 = sensor1_Quater1.Values.data;
sensor1_Quater2 = get(logsout,"sensor1_Quater2");
Sensor1Quater2 = sensor1_Quater2.Values.data;
sensor1_Quater3 = get(logsout,"sensor1_Quater3");
Sensor1Quater3 = sensor1_Quater3.Values.data;
sensor1_Quater4 = get(logsout,"sensor1_Quater4");
Sensor1Quater4 = sensor1_Quater4.Values.data;
sensor1_x = get(logsout,"sensor1_x");
Sensor1X = sensor1_x.Values.data;
sensor1_y = get(logsout,"sensor1_y");
Sensor1Y = sensor1_y.Values.data;
sensor1_z = get(logsout,"sensor1_z");
Sensor1Z = sensor1_z.Values.data;
% concatenating the position and quaterion
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];

flag = get(logsout,"flag");
flagStatus = flag.Values.data;

Motor1_pos = get(logsout,"Motor1_pos_req_qc_mode2");
Motor1Pos = Motor1_pos.Values.data;
Motor2_pos = get(logsout,"Motor2_pos_req_qc_mode2");
Motor2Pos = Motor2_pos.Values.data;
Motor3_pos = get(logsout,"Motor3_pos_req_qc_mode2");
Motor3Pos = Motor3_pos.Values.data;
Motor4_pos = get(logsout,"Motor4_pos_req_qc_mode2");
Motor4Pos = Motor4_pos.Values.data;
uMotorPos = [Motor1Pos';Motor2Pos';Motor3Pos';Motor4Pos'];
%% setup matrix
zu = {};
keptIdxs = 2412:9137;
y_eq = Sensor1(:,9300);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);
u_eq = uMotorPos(:, 9300);
uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorPos(:, keptIdxs) - u_eq; 
%% Multi-stage smoothing for reliable differential calculus
% Parameters
windowSize = 101;  % Must be odd
polynomialOrder = 2;  % Balance between smoothing and feature preservation
medianFilterSize = 7;  % For removing spikes
butterworthOrder = 4;  % Butterworth filter order
cutoffFreq = 0.01;     % Normalized cutoff frequency (0 to 1)
  
for iDim = 1:size(Sensor1,1)
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
load SSM_model.mat
CtrlData5 = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
% pack the trajectory and control input data at the low dimentional as struct
CtrlData5Packed = struct('t', zuSmooth{1,1}', 'u', uCell{1,2}', 'x', CtrlData5{1,2}','y',Yu{1,2}');
%% load new control data 5:trajectory_data_20250611_v2.mat, Gain=1000, spdRegFactor=1
clearvars -except dataDecay CtrlData1Packed CtrlData2Packed CtrlData3Packed CtrlData4Packed CtrlData5Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202506112058.mat

sensor1_Quater1 = get(logsout,"sensor1_Quater1");
time = sensor1_Quater1.Values.Time;
Sensor1Quater1 = sensor1_Quater1.Values.data;
sensor1_Quater2 = get(logsout,"sensor1_Quater2");
Sensor1Quater2 = sensor1_Quater2.Values.data;
sensor1_Quater3 = get(logsout,"sensor1_Quater3");
Sensor1Quater3 = sensor1_Quater3.Values.data;
sensor1_Quater4 = get(logsout,"sensor1_Quater4");
Sensor1Quater4 = sensor1_Quater4.Values.data;
sensor1_x = get(logsout,"sensor1_x");
Sensor1X = sensor1_x.Values.data;
sensor1_y = get(logsout,"sensor1_y");
Sensor1Y = sensor1_y.Values.data;
sensor1_z = get(logsout,"sensor1_z");
Sensor1Z = sensor1_z.Values.data;
% concatenating the position and quaterion
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];

flag = get(logsout,"flag");
flagStatus = flag.Values.data;

Motor1_pos = get(logsout,"Motor1_pos_req_qc_mode2");
Motor1Pos = Motor1_pos.Values.data;
Motor2_pos = get(logsout,"Motor2_pos_req_qc_mode2");
Motor2Pos = Motor2_pos.Values.data;
Motor3_pos = get(logsout,"Motor3_pos_req_qc_mode2");
Motor3Pos = Motor3_pos.Values.data;
Motor4_pos = get(logsout,"Motor4_pos_req_qc_mode2");
Motor4Pos = Motor4_pos.Values.data;
uMotorPos = [Motor1Pos';Motor2Pos';Motor3Pos';Motor4Pos'];
%% setup matrix
zu = {};
keptIdxs = 3330:13311;
y_eq = Sensor1(:,13566);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);
u_eq = uMotorPos(:, 13566);
uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorPos(:, keptIdxs) - u_eq; 
%% Multi-stage smoothing for reliable differential calculus
% Parameters
windowSize = 101;  % Must be odd
polynomialOrder = 2;  % Balance between smoothing and feature preservation
medianFilterSize = 7;  % For removing spikes
butterworthOrder = 4;  % Butterworth filter order
cutoffFreq = 0.01;     % Normalized cutoff frequency (0 to 1)
  
for iDim = 1:size(Sensor1,1)
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
load SSM_model.mat
CtrlData6 = projectTrajectories(ssm_model.IMInfoCtrl, Yu); % Project down
% pack the trajectory and control input data at the low dimentional as struct
CtrlData6Packed = struct('t', zuSmooth{1,1}', 'u', uCell{1,2}', 'x', CtrlData6{1,2}','y',Yu{1,2}');
%% save data
lowDimCtrl{1} = CtrlData1Packed;
lowDimCtrl{2} = CtrlData2Packed;
lowDimCtrl{3} = CtrlData3Packed;
lowDimCtrl{4} = CtrlData4Packed;
lowDimCtrl{5} = CtrlData5Packed;
lowDimCtrl{6} = CtrlData6Packed;
% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('processedData', sprintf('capacityKoopmanObs_%s.mat', current_date));

% Save the trajectory data
save(fileName, 'lowDimCtrl');
%%
clearvars -except lowDimCtrl dataDecay;
clc
close all
%% Set test and train trajectories out of all the trajectories
indOut = [];
indTest = [8 6];
indTrain = setdiff(1:length(dataDecay),union(indTest,indOut));
decayNuminTrain = length(indTrain);
decayNuminTest = length(indTest);
%% reshape a data file with name like: softrobot_train-13_val-4.mat
for i = 1:length(dataDecay)
    dataDecay{1,i}.u = zeros(size(dataDecay{1,i}.t, 1),4);
    dataDecayResmp{1,i} = resample(dataDecay{1,i},0.01);
end
train = dataDecayResmp(:,indTrain);
val = dataDecayResmp(:,indTest);
%% Learn control matrix
yuxData = cell(size(lowDimCtrl,1),3);
for i = 1:size(lowDimCtrl,1)
    yuxData{i,1} = lowDimCtrl{i,1};
    yuxData{i,2} = lowDimCtrl{i,2};
    yuxData{i,3} = lowDimCtrl{i,3};
    yuxData{i,4} = lowDimCtrl{i,4};
    yuxData{i,5} = lowDimCtrl{i,5};
    yuxData{i,6} = lowDimCtrl{i,6};
end
%% Mix the control data
% Get half the length of each segment
halfLength1 = floor(size(yuxData{1,1}.t, 1) / 2);
halfLength2 = floor(size(yuxData{1,2}.t, 1) / 2);

% Extract half of the first segment
x1_half1 = yuxData{1,1}.x(1:halfLength1,:);
u1_half1 = yuxData{1,1}.u(1:halfLength1,:);
y1_half1 = yuxData{1,1}.y(1:halfLength1,:);
x1_half2 = yuxData{1,1}.x(halfLength1+1:end,:);
u1_half2 = yuxData{1,1}.u(halfLength1+1:end,:);
y1_half2 = yuxData{1,1}.y(halfLength1+1:end,:);

% Extract half of the second segment
x2_half1 = yuxData{1,2}.x(1:halfLength2,:);
u2_half1 = yuxData{1,2}.u(1:halfLength2,:);
y2_half1 = yuxData{1,2}.y(1:halfLength2,:);
x2_half2 = yuxData{1,2}.x(halfLength1+1:end,:);
u2_half2 = yuxData{1,2}.u(halfLength1+1:end,:);
y2_half2 = yuxData{1,2}.y(halfLength1+1:end,:);

% Mix and concatenate
mixedData1 = struct('u', [u1_half1; u2_half1], 'y', [y1_half1; y2_half1],'x',[x1_half1; x2_half1]);
mixedData1.t = (0.01:0.01:size(mixedData1.u,1)*0.01)';

mixedData2 = struct('u', [u1_half2; u2_half2], 'y', [y1_half2; y2_half2],'x',[x1_half2; x2_half2]);
mixedData2.t = (0.01:0.01:size(mixedData2.u,1)*0.01)';

mixedData3 = lowDimCtrl{1, 4};
mixedData4 = lowDimCtrl{1, 5};
mixedData5 = lowDimCtrl{1, 6};

dataCtrl{1} = mixedData1;
dataCtrl{2} = mixedData2;
dataCtrl{3} = mixedData3;
dataCtrl{4} = mixedData4;
dataCtrl{5} = mixedData5;
%% Set test and train trajectories out of the controlled trajectories
indTest = 5;
indTrain = setdiff(1:length(dataCtrl),indTest);
%% 
train = [train,dataCtrl(indTrain)];
val = [val,dataCtrl(indTest)];
% Get the lengths of train and val cell arrays
a = length(train);
b = length(val);

% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('processedData', sprintf('softKpmObs_capacity_train-%d_val-%d_%s.mat', a, b, current_date));

% Save the cell arrays to the .mat file
save(fileName, 'train', 'val', 'decayNuminTrain', "decayNuminTest");