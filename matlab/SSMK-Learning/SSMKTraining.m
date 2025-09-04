% SoftMnp_train
%
% Creates a sysid class and walks through all of the steps of training a
% koopman-based model from data, validating its performance, and saving it (if desired)
% Last update: June/11/2025 by Yongjun Yan
% Revision history: 


%% gather training data (need to prepare data file before running this) (load softrobot_train-15_val-3.mat)
clc
clear all
close all
% load in data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'processedData/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );


%% construct sysid class
ksysid = ksysid( data4sysid, ...
        'model_type' , 'linear' ,...    % model type (linear or nonlinear)
        'obs_type' , { 'gaussian' } ,...    % type of basis functions:poly, fourier, fourier-sparser, gaussian, hermite
        'obs_degree' , [ 2 ] ,...       % "degree" of basis functions
        'snapshots' , Inf ,...          % Number of snapshot pairs
        'lasso' , [ Inf ] ,...           % L1 regularization term
        'delays' , 1 );                 % Numer of state/input delays

    
%% train model(s)
ksysid = ksysid.train_models;

%% validate model(s)
% could also manually do this for one model at a time

results = cell( size(ksysid.candidates) );    % store results in a cell array
err = cell( size(ksysid.candidates) );    % store error in a cell array 

if iscell(ksysid.candidates)
    for i = 1 : length(ksysid.candidates)
        [ results{i} , err{i} ] = ksysid.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = ksysid.valNplot_model;
end
load processedData/figure8SSMKOBsEmbd_2025-08-07.mat
yDataPredHoriz = ksysid.valOverHorizion_lowDim(ObsEmdbCtrl{1,3});
%% check at the observable space
load SSM_model.mat
for i = 1:size(yDataPredHoriz,1)
    yDataPredHoriz{i,6} = liftTrajectories(ssm_model.IMInfoCtrl, {yDataPredHoriz{i,3}', yDataPredHoriz{i,1}'});
    outdofsDelay = [2 3];
    trajectory_est = yDataPredHoriz{i,6}{1,2}(outdofsDelay,:)';
    trajectory_ref = yDataPredHoriz{i,4}(2:end,outdofsDelay);
    [rmse_pos,  pos_errors,] = computeTrajectoryError(trajectory_est, trajectory_ref);
    % Initialize error structure
    err = struct('rmsePos', [], 'rmseRot', [], 'posErrorsc', [], 'rotErrorsc', []);
    err.rmsePos = rmse_pos;
    err.posErrors = pos_errors;
    yDataPredHoriz{i,7} = err;
end
%%
[ yDataPredHoriz ] = valOverHorizion_plot(yDataPredHoriz, ObsEmdbCtrl{1,3}', ksysid.valdata{1,1});
%%
% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('results', sprintf('dataCapaPredHorizSSMKoopman_%s.mat', current_date));

% Save the trajectory data
save(fileName, 'yDataPredHoriz');
%%
function [rmse_pos,  pos_errors] = computeTrajectoryError(trajectory_est, trajectory_ref)
    % trajectory_est, trajectory_ref: Nx7 matrices [x y z qx qy qz qw]

    % Ensure input sizes match
    assert(all(size(trajectory_est) == size(trajectory_ref)), 'Trajectory sizes must match');
    
    N = size(trajectory_est, 1);
    pos_errors = zeros(N, 1);
    rot_errors = zeros(N, 1);
    diffs = zeros(N, 1);

    for i = 1:N
        % Extract positions
        p_est = trajectory_est(i,:);
        p_ref = trajectory_ref(i,:);

        % Position error (Euclidean)
        pos_errors(i) = norm(p_est - p_ref);
    end

    % Compute RMSE
    rmse_pos = sqrt(mean(pos_errors.^2));
    
    fprintf('Position RMSE: %.4f mm\n', rmse_pos);
end

% Validate the derived model over short predictive horizion
function [ yDataPredHoriz ] = valOverHorizion_plot(yDataPredHoriz, validateData, valData)
     % validation at the reduced coordinates
    idxPlotCL = 50;
    % Plot time series
    customFigure('subPlot',[2 1]);
    subplot(2,1,1); title ('Reduced dynamics')
    for i = 2:3
        subplot(2,1,i-1);
        xlabel('t','Interpreter','latex');
        if i == 2
            ylabel('Y [mm]');
        else
            ylabel('Z [mm]');
        end
    end
    % for iPlt = 1:7
    %     subplot(7,1,iPlt);
    %     xlabel('$t$','Interpreter','latex');
    %     ylabel(labelsRD(iPlt),'Interpreter','latex');
    % end
    for iObs = 2:3
        subplot(2,1,iObs-1);
        plot(valData.t, validateData(iObs,1:size(valData.t,1)),'Linewidth',1);
    end
    for iObs = 2:3
        subplot(2,1,iObs-1);
        for iTraj=1:idxPlotCL:size(yDataPredHoriz,1)
            plot(yDataPredHoriz{iTraj,6}{1,1}, yDataPredHoriz{iTraj,6}{1,2}(iObs,:),'Linewidth',2);
        end
    end
end