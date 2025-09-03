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
%% validate over short predictive horizon
yDataPredHoriz = ksysid.valOverHorizion_obs;
%% check at the observable space
load SSM_model.mat
for i = 1:size(yDataPredHoriz,1)
    trajectory_est = yDataPredHoriz{i,1};
    trajectory_ref = yDataPredHoriz{i,4}(2:end,:);
    [rmse_pos, rmse_rot, pos_errors, rot_errors] = computeTrajectoryError(trajectory_est, trajectory_ref);
    % Initialize error structure
    err = struct('rmsePos', [], 'rmseRot', [], 'posErrorsc', [], 'rotErrorsc', []);
    err.rmsePos = rmse_pos;
    err.rmseRot = rmse_rot;
    err.posErrors = pos_errors;
    err.rotErrors = rot_errors;
    yDataPredHoriz{i,7} = err;
end
[ yDataPredHoriz ] = valOverHorizion_plot(yDataPredHoriz, data4sysid.val{1,3});
%%
% Get the current date
current_date = datetime('now','TimeZone','local','Format','yyyy-MM-dd'); % Format: YYYYMMDD

% Create the dynamic file name with the folder path
fileName = fullfile('results', sprintf('dataCapaPredHorizKoopmanObs_%s.mat', current_date));

% Save the trajectory data
save(fileName, 'yDataPredHoriz');
%%
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

% Validate the derived model over short predictive horizion
function [ yDataPredHoriz ] = valOverHorizion_plot(yDataPredHoriz, validateData)
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
        plot(validateData.t, validateData.y(:,iObs),'Linewidth',1);
    end
    for iObs = 1:7
        subplot(7,1,iObs);
        for iTraj=1:idxPlotCL:size(yDataPredHoriz,1)
            plot(yDataPredHoriz{iTraj,3}, yDataPredHoriz{iTraj,1}(:,iObs),'Linewidth',2);
        end
    end
end