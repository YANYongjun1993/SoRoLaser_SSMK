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
%% Delay-embedding
SSMDim = 4;
overEmbed = 32;
[yData, opts_embd] = coordinatesEmbeddingControl(oDataTrunc, SSMDim, 'OverEmbedding', overEmbed);
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
    IMInfo.chart.map = @(x) transpose(Vde)*x;
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
IMInfo.chart.map = @(x) transpose(Vde)*x;
IMInfo = struct('chart', IMInfo.chart, 'parametrization', IMInfo.parametrization);

for i = 1:nPoints
    % Assign struct to the i-th cell
    dataDecay{i} = struct('t', zDataDelay{i,1}', 'y', zDataDelay{i,2}', 'x', zDataDelay{i,2}');
end
%% Export mapping coefficients and input matrix
ssm_model = struct();
ssm_model.sizeX = 4;
ssm_model.sizeU = 4;
ssm_model.SSMDim = SSMDim;
ssm_model.IMInfo = IMInfo;
save('SSM_model.mat', 'ssm_model', '-v7');
% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('processedData', sprintf('koopmanDecayObs_%s.mat', current_date));

% Save the trajectory data
save(fileName, 'dataDecay');
%% load new control data 1: trajectory_data_20250801_v2.mat, spdRegFactor = 2; gain = 800
clearvars -except dataDecay;
clc
close all
% Load new controller data with random inputs
load rawData/202508011703.mat

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
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];

yMotor_cmd = get(logsout,"yMotor");
yMotor = yMotor_cmd.Values.data;
zMotor_cmd = get(logsout,"zMotor");
zMotor = zMotor_cmd.Values.data;
uMotorCmd = [yMotor';zMotor'];
%% setup matrix
zu = {};
keptIdxs = 6415:18416;
y_eq = Sensor1(:,6272);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);
uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorCmd(:, keptIdxs); 

% pack the trajectory and control input data at the low dimentional as struct
CtrlData1Packed = struct('t', zu{1,1}', 'u', uCell{1,2}', 'y', zu{1,2}','x',zu{1,2}');


%% load new control data 2:  trajectory_data_20250801_v2.mat, spdRegFactor = 2; gain = 1100
clearvars -except dataDecay ObsEmdbCtrl CtrlData1Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202508011714.mat
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
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];

flag = get(logsout,"flag");
flagStatus = flag.Values.data;

yMotor_cmd = get(logsout,"yMotor");
yMotor = yMotor_cmd.Values.data;
zMotor_cmd = get(logsout,"zMotor");
zMotor = zMotor_cmd.Values.data;
uMotorCmd = [yMotor';zMotor'];
%% setup matrix
zu = {};
keptIdxs = 19575:31572;
y_eq = Sensor1(:,19418);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);

uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorCmd(:, keptIdxs); 
% pack the trajectory and control input data at the low dimentional as struct
CtrlData2Packed = struct('t', zu{1,1}', 'u', uCell{1,2}', 'y', zu{1,2}','x',zu{1,2}');

%% load new control data 3: trajectory_data_20250801_v1.mat, gain=1100; spdRegFactor = 2
clearvars -except dataDecay ObsEmdbCtrl CtrlData1Packed CtrlData2Packed;
clc
close all
% Load new controller data with random inputs
load rawData/202508011726.mat
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
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];

flag = get(logsout,"flag");
flagStatus = flag.Values.data;

yMotor_cmd = get(logsout,"yMotor");
yMotor = yMotor_cmd.Values.data;
zMotor_cmd = get(logsout,"zMotor");
zMotor = zMotor_cmd.Values.data;
uMotorCmd = [yMotor';zMotor'];
%% setup matrix
zu = {};
keptIdxs = 8775:20772;
y_eq = Sensor1(:,8467);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);

uCell{1,1} = zu{1,1}; uCell{1,2} = uMotorCmd(:, keptIdxs); 
% pack the trajectory and control input data at the low dimentional as struct
CtrlData3Packed = struct('t', zu{1,1}', 'u', uCell{1,2}', 'y', zu{1,2}','x',zu{1,2}');

lowDimCtrl{1} = CtrlData1Packed;
lowDimCtrl{2} = CtrlData2Packed;
lowDimCtrl{3} = CtrlData3Packed;
%%
clearvars -except lowDimCtrl dataDecay;
clc
close all

%% Learn control matrix
% Set test and train trajectories out of the controlled trajectories
indTest = 3;
indTrain = setdiff(1:length(lowDimCtrl),indTest);
%% 
train = [lowDimCtrl(indTrain)];
val = [lowDimCtrl(indTest)];
% Get the lengths of train and val cell arrays
a = length(train);
b = length(val);

% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('processedData', sprintf('softKoopman_figure8_train-%d_val-%d_%s.mat', a, b, current_date));

% Save the cell arrays to the .mat file
save(fileName, 'train', 'val');