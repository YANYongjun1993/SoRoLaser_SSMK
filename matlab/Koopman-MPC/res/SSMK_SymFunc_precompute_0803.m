clc
clear all
% load in data file(s) n-4_m-4_del-1_2025-07-05_16-56.mat
[ datafile_name , datafile_path ] = uigetfile( 'E:\Research\SRIL\14th\TroRevision\HardwareExperiemnt\SSMKTraining\systems\fromData\*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );
% load parameterized SSM model SSM_model.m
load E:\Research\SRIL\14th\TroRevision\HardwareExperiemnt\SSMKTraining\SSM_model.mat

[f]=generate_functions(data4sysid.ksysid, ssm_model.IMInfo);

function [f ]=generate_functions( data4sysid,  ssm_model)
%%-------------------------------------------------------------------------        
% Get users' functions for generating sub-functions of users' optimization 
% problems from the trained koopman-model in .\Koopman_Model_Learning folder.
%
% Last update: July/21/2025 by Yongjun Yan
% Revision history: generate the symbolic expression for the Koopman-based
% MPC algorithm
%  x     : The state at the lifted space.
%  u     : The control variables.
%  y     : The state at the low-dimensional space.
%  z    : The delay-embeded observables.
%  zeta : The low dimensional state with delays.

syms x u y z zeta real
for n = 1:1:data4sysid.params.N
    syms (sprintf('x%01d',n));
    x(n) = sym(sprintf('x%01d',n),'real');
end
for m = 1:1:data4sysid.params.m
    syms (sprintf('u%01d',m));
    u(m) = sym(sprintf('u%01d',m),'real');
end
for o = 1:1:data4sysid.params.n
    syms (sprintf('y%01d',o));
    y(o) = sym(sprintf('y%01d',o),'real');
end
for i = 1:1:data4sysid.params.nzeta
    syms (sprintf('zeta%01d',i));
    zeta(i) = sym(sprintf('zeta%01d',i),'real');
end
for i = 1:1:size(ssm_model.parametrization.tangentSpaceAtOrigin,1)
    syms (sprintf('z%01d',i));
    z(i) = sym(sprintf('z%01d',i),'real');
end

scaleup_u = data4sysid.scaleup.u(u);
scaleup_y = data4sysid.scaleup.y(y);
scaledown_u = data4sysid.scaledown.u(u);
scaledown_y = data4sysid.scaledown.y(y);
lift_full = data4sysid.lift.full(zeta');
f = data4sysid.model.A*x' + data4sysid.model.B*u';
D_lowDim = data4sysid.model.C*(x');
% chart map
chart = ssm_model.chart.map(z'); % Project down
% parameterization map
param = ssm_model.parametrization.map(y');

if isempty(scaleup_u)
    scaleup_u = zeros(0,1);
    matlabFunction(scaleup_u, 'file', 'scaleup_u','vars',[u]);
else
    matlabFunction(scaleup_u, 'file', 'scaleup_u','vars',{u}); 
end

if isempty(scaleup_y)
    scaleup_y = zeros(0,1);
    matlabFunction(scaleup_y, 'file', 'scaleup_y','vars',[y]);
else
    matlabFunction(scaleup_y, 'file', 'scaleup_y','vars',{y}); 
end

if isempty(scaledown_u)
    scaledown_u = zeros(0,1);
    matlabFunction(scaledown_u, 'file', 'scaledown_u','vars',[u]);
else
    matlabFunction(scaledown_u, 'file', 'scaledown_u','vars',{u}); 
end

if isempty(f)
    f = zeros(0,1);
    matlabFunction(f, 'file', 'valf','vars',[x, u]);
else
    matlabFunction(f, 'file', 'valf','vars',{x, u}); 
end

if isempty(scaledown_y)
    scaledown_y = zeros(0,1);
    matlabFunction(scaledown_y, 'file', 'scaledown_y','vars',[y]);
else
    matlabFunction(scaledown_y, 'file', 'scaledown_y','vars',{y}); 
end

if isempty(lift_full)
    lift_full = zeros(0,1);
    matlabFunction(lift_full, 'file', 'lift_full','vars',[zeta]);
else
    matlabFunction(lift_full, 'file', 'lift_full','vars',{zeta}); 
end

if isempty(D_lowDim)
    D_lowDim = zeros(0,1);
    matlabFunction(D_lowDim, 'file', 'D_lowDim','vars',[x]);
else
    matlabFunction(D_lowDim, 'file', 'D_lowDim','vars',{x}); 
end

if isempty(chart)
    chart = zeros(0,1);
    matlabFunction(chart, 'file', 'chart','vars',[z]);
else
    matlabFunction(chart, 'file', 'chart','vars',{z}); 
end

if isempty(param)
    param = zeros(0,1);
    matlabFunction(param, 'file', 'param','vars',[y]);
else
    matlabFunction(param, 'file', 'param','vars',{y}); 
end

return
end