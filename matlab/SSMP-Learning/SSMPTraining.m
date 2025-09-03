clearvars;
clc
clear all
close all
%% load data
load processedData/SSMPDecayObs_2025-08-13.mat
load processedData/figure8SSMP_2025-08-13.mat
for i = 1:length(lowDimCtrl)
    lowDimCtrl{i} = rmfield(lowDimCtrl{i},'x');
end
load SSM_model.mat
%% Set test and train trajectories out of all the trajectories
indOut = [3 9 13 15];
indTest = [2 8 11];
indTrain = setdiff(1:length(dataDecay),union(indTest,indOut));
%% scale all the date to [-1,1]
zDataMerged = merge_trials(dataDecay,lowDimCtrl);
[zDataDecayNorm, zDataCtrlNorm, obj] = get_scale(zDataMerged, dataDecay, lowDimCtrl, ssm_model);
%% save data
% dataNorm{1} = zDataDecayNorm;
% dataNorm{2} = zDataCtrlNorm;
% % Get the current date
% current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD
% 
% % Create the dynamic file name with the folder path
% fileName = fullfile('processedData', sprintf('dataNormPoly_%s.mat', current_date));
% 
% % Save the trajectory data
% save(fileName, 'dataNorm');
%% Calibrate reduced dynamics in the new reduced coordinates
ROMOrderObs = 1;
Nregvals = 9;
RDInfoCtrl = IMDynamicsFlowCV(zDataDecayNorm(indTrain,:), 'R_PolyOrd', ROMOrderObs, 'l_vals', 0.001,'style', 'default');
% RDInfo = IMDynamicsFlowCV(yDataTrunc, 'R_PolyOrd', ROMOrder, 'n_folds', Nfolds, 'l_vals', logspace(-6,0,Nregvals),'style', 'default');
RDInfoCtrl.eigenvaluesLinPartFlow(1:ssm_model.SSMDim)
Rauton = RDInfoCtrl.reducedDynamics.map;
reDyn = Rauton(zDataDecayNorm{1,2});
disp('Errors for the dynamics on reduced-coordinates:')
zRec = advectRD(RDInfoCtrl, zDataDecayNorm);
normedTrajDist = computeTrajectoryErrors(zRec, zDataDecayNorm);
meanErrorDynTrain = mean(normedTrajDist(indTrain))*100
meanErrorDynTest = mean(normedTrajDist(indTest))*100
yDataRec = liftTrajectories(ssm_model.IMInfoCtrl, zRec);

% Plot settings
vPlot = 1:length(indTest);
nDims = ssm_model.SSMDim;  % Number of dimensions to plot
customFigure('subPlot',[nDims 1]);

% Create subplots dynamically
for i = 1:nDims
    subplot(nDims,1,i);
    colororder(cool(length(vPlot)));
    xlabel('$t$','Interpreter','latex');
    ylabel(sprintf('$X_%d$',i),'Interpreter','latex');
end
% Plot trajectories for the reconstrunction accuracy comparision at the
% reduced coordinates
j = 1;
for iTraj = vPlot
    % Check for NaN values in SSMDim dimensions
    if ~any(isnan(zRec{iTraj,2}(1:nDims,:)), 'all')
        % Plot each dimension
        for i = 1:nDims
            subplot(nDims,1,i);
            % Plot reconstructed trajectory
            plot(zRec{iTraj,1}, zRec{iTraj,2}(i,:), 'k:', 'LineWidth', 1);
            hold on;
            % Plot original trajectory
            plot(zDataDecayNorm{iTraj,1}, zDataDecayNorm{iTraj,2}(i,:), 'LineWidth', 1);
            hold on
        end
    else
        bad(j) = iTraj;
        j = j + 1;
    end
end
%% 
indTest = [3];
indTrain = setdiff(1:size(lowDimCtrl,2),indTest);
%% Fit control matrices
[B,Br,regErrorsNorm,regErrorsAvg] = fitControlMatricesCV(zDataCtrlNorm(indTrain,:),RDInfoCtrl,ssm_model.IMInfoCtrl,0.0)
% Update models
Rauton = RDInfoCtrl.reducedDynamics.map;
Monom = RDInfoCtrl.reducedDynamics.phi;
R = @(t,x,u) Rauton(x) + B*u + Br*kron(Monom(x),u);
Rctrl = @(x,u) Rauton(x) + B*u + Br*kron(Monom(x),u);
%% validation of the whole trajectory at the reduced coordinates
% check dynamics
ti = zDataCtrlNorm{indTest(1),1};
xi = zDataCtrlNorm{indTest(1),2};
ui = zDataCtrlNorm{indTest(1),3};
uFun = @(t) transpose(interp1(ti, transpose(ui), t, 'pchip', 'extrap'));
% opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-2);
[tiPred,xiPred] = ode45(@(t,x) R(t,x,uFun(t)),ti,xi(:,1)); 
tiPred = transpose(tiPred); xiPred = transpose(xiPred);
labelsRD = ["$\alpha_1$", "$\alpha_2$", "$\alpha_3$", "$\alpha_4$"];
% Plot time series
customFigure('subPlot',[2 1]);
subplot(2,1,1); title ('Reduced dynamics')
for iPlt = 1:2
    subplot(2,1,iPlt);
    xlabel('$t$','Interpreter','latex');
    ylabel(labelsRD(iPlt),'Interpreter','latex');
end
for iObs = 1:2
    subplot(2,1,iObs);
    plot(ti, xi(iObs, :),'Linewidth',1);
end
for iObs = 1:2
    subplot(2,1,iObs);
    plot(tiPred, xiPred(iObs, :),'Linewidth',1);
end

%% Validate short time Control model
load processedData/figure8SSMP_2025-08-13.mat
load processedData/figure8SSMPObs_2025-08-13.mat
xObswo = ObsEmdbCtrl{1,indTest(1)}(1:size(ui,2),:)';
% xObswo = lowDimCtrl{1,indTest(1)}.x';
Nhoriz = 11;
nEval = length(ti)-Nhoriz;
yDataPredHoriz = cell(1,5);
iRow = 1;
for iTime = 1:nEval
    uHoriz = ui(:,iTime:iTime+Nhoriz);
    tHoriz = ti(iTime:iTime+Nhoriz);
    xHoriz = xObswo(:,iTime:iTime+Nhoriz);
    uFun = @(t) transpose(interp1(tHoriz, transpose(uHoriz), t));
    [iTimeRec, ietaRec] = ode15s(@(t,x) R(t,x,uFun(t)),tHoriz,xi(:,iTime));
    % Lift
    ietaObs = obj.scaleup.y(ietaRec);
    iyRec = liftTrajectories(ssm_model.IMInfoCtrl, {iTimeRec, transpose(ietaObs)});
    
    yDataPredHoriz{iRow,1} = ti(iTime:iTime+Nhoriz);
    yDataPredHoriz{iRow,2} = ietaRec;
    yDataPredHoriz{iRow,3} = ietaObs;
    yDataPredHoriz{iRow,4} = iyRec{1,2};
    yDataPredHoriz{iRow,5} = xHoriz;
    iRow = iRow + 1;
end
%% validation at the reduced coordinates
idxPlotCL = 50;
labelsRD = ["$\alpha_1$", "$\alpha_2$", "$\alpha_3$", "$\alpha_4$"];
% Plot time series
customFigure('subPlot',[2 1]);
subplot(2,1,1); title ('Reduced dynamics')
for iPlt = 1:2
    subplot(2,1,iPlt);
    xlabel('$t$','Interpreter','latex');
    ylabel(labelsRD(iPlt),'Interpreter','latex');
end
for iObs = 1:2
    subplot(2,1,iObs);
    plot(ti, xi(iObs, :),'Linewidth',1);
end
for iObs = 1:2
    subplot(2,1,iObs);
    for iTraj=1:idxPlotCL:size(yDataPredHoriz,1)
        plot(yDataPredHoriz{iTraj,1}, yDataPredHoriz{iTraj,2}(:,iObs),'Linewidth',2);
    end
end
%%
for i = 1:size(yDataPredHoriz,1)
    trajectory_est = yDataPredHoriz{i,4}';
    trajectory_ref = yDataPredHoriz{i,5}';
    [rmse_pos, rmse_rot, pos_errors, rot_errors] = computeTrajectoryError(trajectory_est, trajectory_ref);
    % Initialize error structure
    err = struct('rmsePos', [], 'rmseRot', [], 'posErrorsc', [], 'rotErrorsc', []);
    err.rmsePos = rmse_pos;
    err.rmseRot = rmse_rot;
    err.posErrors = pos_errors;
    err.rotErrors = rot_errors;
    yDataPredHoriz{i,6} = err;
end
%%
[ yDataPredHoriz ] = valOverHorizion_plot(yDataPredHoriz, lowDimCtrl{1,indTest(1)}, ObsEmdbCtrl{1,indTest(1)}(1:size(ui,2),:));
%%
% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('results', sprintf('dataPredHorizSSMPoly_%s.mat', current_date));

% Save the trajectory data
save(fileName, 'yDataPredHoriz');
%% Export mapping coefficients and input matrix
% Continuous-time model

% TODO: Save either known input matrix or learned input matrix

Ts = 0.01;

py_data = struct();

ssm_model.Br = Br;
ssm_model.B = B;
ssm_model.R= Rctrl;
ssm_model.RDInfoCtrl = RDInfoCtrl;

ssm_model.Ts = Ts;
py_data.model = ssm_model;

params = struct();
params.SSM_order = 1;
params.ROM_order = 1;
params.state_dim = 2;
params.input_dim = 2;
params.output_dim = 7;
params.delays = 32;
py_data.params = params;
py_data.scale = obj;

save('SSMP_model_0814.mat', 'py_data', '-v7')

%% Validate the derived model over short predictive horizion
function [ yDataPredHoriz ] = valOverHorizion_plot(yDataPredHoriz, validateDataLowDim, validateDataObs)
     % validation at the reduced coordinates
    idxPlotCL = 50;
    labelsRD = ["$\alpha_1$", "$\alpha_2$", "$\alpha_3$", "$\alpha_4$"];
    % Plot time series
    customFigure('subPlot',[7 1]);
    xyz = 'xyz'; % Create string for axis labels
    subplot(7,1,1); title ('Reduced dynamics')
    for i = 1:7
        subplot(7,1,i);
        xlabel('t','Interpreter','latex');
        if i <= 3
            ylabel(sprintf('%s-axis [mm]', xyz(i)),'Interpreter','latex');
        else
            ylabel(sprintf('q%d',i-3),'Interpreter','latex');
        end
    end
    % for iPlt = 1:7
    %     subplot(7,1,iPlt);
    %     xlabel('$t$','Interpreter','latex');
    %     ylabel(labelsRD(iPlt),'Interpreter','latex');
    % end
    for iObs = 1:7
        subplot(7,1,iObs);
        plot(validateDataLowDim.t, validateDataObs(:,iObs),'Linewidth',1);
    end
    for iObs = 1:7
        subplot(7,1,iObs);
        for iTraj=1:idxPlotCL:size(yDataPredHoriz,1)
            plot(yDataPredHoriz{iTraj,1}, yDataPredHoriz{iTraj,4}(iObs,:),'Linewidth',2);
        end
    end
end
%% compute Trajectory Error: position and orientation
function [rmse_pos, rmse_angle_deg, pos_errors, rot_errors] = computeTrajectoryError(trajectory_est, trajectory_ref)
    % trajectory_est, trajectory_ref: Nx7 matrices [x y z qx qy qz qw]

    % Ensure input sizes match
    assert(all(size(trajectory_est) == size(trajectory_ref)), 'Trajectory sizes must match');
    
    N = size(trajectory_est, 1);
    pos_errors = zeros(N, 1);
    rot_errors = zeros(N, 1);
    diffs = zeros(N, 1);

    for i = 1:N
        % Extract positions
        p_est = trajectory_est(i, 1:3);
        p_ref = trajectory_ref(i, 1:3);

        % Position error (Euclidean)
        pos_errors(i) = norm(p_est - p_ref);

        % Extract orientations (quaternions)
        q_est = trajectory_est(i, 4:7);
        q_ref = trajectory_ref(i, 4:7);

        % Quaternion to rotation matrix
        R_est = quat2rotm(q_est);
        R_ref = quat2rotm(q_ref);

        % Extract Z-axis direction vector (third column)
        z_est = R_est(:, 3);
        z_ref = R_ref(:, 3);

        % Angle (angular error)
        cos_theta = dot(z_est, z_ref) / (norm(z_est) * norm(z_ref));
        cos_theta = min(1, max(-1, cos_theta));  % Avoid numerical errors from exceeding valid range
        theta = acos(cos_theta);  % radian
        rot_errors(i) = theta;

        % Euclidean distance error
        diffs(i) = norm(z_est - z_ref);
    end

    % Compute RMSE
    rmse_pos = sqrt(mean(pos_errors.^2));
    rmse_angle_deg = rad2deg(sqrt(mean(rot_errors.^2)));
    rmse_vector = sqrt(mean(diffs.^2));
    
    fprintf('Position RMSE: %.4f mm\n', rmse_pos);
    fprintf('Rotation RMSE: %.4f °\n', rmse_angle_deg);
    fprintf('Z-axis vector error RMSE: %.4f radian\n', rmse_vector);
end