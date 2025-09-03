clc
clear all
% load in data file(s) n-4_m-4_del-1_2025-07-05_16-56.mat
[ datafile_name , datafile_path ] = uigetfile( 'E:\Research\SRIL\14th\TroRevision\HardwareExperiemnt\Koopman\KoopmanTraining\systems\fromData\*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );
MPC.horizon = 10;
var.costRun = 0.1*diag([0 1 1 0 0 0 0]);
var.costTer = 100*diag([0 1 1 0 0 0 0]);
var.costInp = 0.01;
var.Inb = [-1100, 1100];
var.InSlp = 0.5e-2;
var.InSmt = 2e-1;
% params.nd = 1;
% params.m = 4;
% params.n = 4;
% isolate the control variables from the delay-embedded observables
% refIso = zeros(2,size(paramMapSSM,1));
% indiRef = [size(paramMapSSM,1)-5, size(paramMapSSM,1)-4];
% refIso(1,indiRef(1)) = 1;
% refIso(2,indiRef(2)) = 1;
[H, G, D, c, L, M]=generate_functions(MPC,data4sysid, var);
save ('MtxMPC.mat','H','G','D','c','L','M')


function [H, G, D, c, L, M]=generate_functions(MPC, data4sysid, var)
%%-------------------------------------------------------------------------        
% Get users' functions for generating sub-functions of users' optimization 
% problems from the trained koopman-model in E:\Research\SRIL\14th\
% TroRevision\HardwareExperiemnt\SSMKTraining\systems\fromData folder.
%
% Last update: July/21/2025 by Yongjun Yan
% Revision history: generate the symbolic expression for the Koopman-based
% MPC algorithm
%% ------------------------------------------------------------------------
% Do not modify below this line 
%% ------------------------------------------------------------------------
% Users' functions: Please provide your functions
% Users should provide following 5 functions in the symbolic form.
% Phi:    terminal cost function
% L:      incremental cost function
% f:      equations of system dynamics 
% C:      input and state mixed constraints
% Cbar:   state only constraints
% Example:
% if your cost function 
%  J=0.5*(y(N)-y_r(N))'*Q_N*(y(N)-y_r(N)) + 0.5*\Sum_(k=0)^(k=N-1)((y(N)-y_r(N))'*Q*(y(N)-y_r(N)) + u(k)'*R*u(k),
% and your eqautaion of system dynamics:
%             h(k+1) = obj.model.A * h(k) + obj.model.B * u(k);
%             z(k) = obj.model.C * h(k);
% and mapping from low-dimensional space to the delay-embeded observables:
%             \tilde{y(k)} = paramMapSSM * z(k);
% and isolate the control variables from the delay-embeded observables for
%             cost function formulation:
%             y(k) = refIso * \tilde{y(k)};
% and constraints for control variables:
%             u(k) \in var.Inb.
% -------------------------------------------------------------------------
% Pre-calculate the matrix 'H','G','D','c','L','M', to cast the MPC problem
% into a standard QP probelem with Cost function:
%            U'HU + ( h0'G + Yr'D )U
% and constrints:
%            LU + M h_0 <= c,
%            h0 = Koopman(\tilde{y(k)}).
    
model = data4sysid.ksysid.model;

% A
N = size(model.A,1);
A = sparse( N*(MPC.horizon+1) , N );
for i = 0 : MPC.horizon
    A( (N*i + 1) : N*(i+1) , : ) = model.A^i ;
end

% B
Bheight = N*(MPC.horizon+1);
Bcolwidth = size(model.B,2);
Bcol = sparse( Bheight , Bcolwidth );    % first column of B matrix
for i = 1 : MPC.horizon
    Bcol( (N*i + 1) : N*(i+1) , : ) = model.A^(i-1) * model.B ;
end

Lshift = spdiags( ones( N*MPC.horizon , 1 ) , -N , N*(MPC.horizon+1) , N*(MPC.horizon+1) );    % lower shift operator

Bwidth = size(model.B,2)*(MPC.horizon);    % total columns in B matrix
Bblkwidth = MPC.horizon;   % how many Bcol blocks wide B matrix is
B = spalloc( Bheight , Bwidth , floor(Bheight * Bwidth / 2) ); % initialze sparse B matrix
B(: , 1:Bcolwidth) = Bcol;
for i = 2 : Bblkwidth
    B(: , (i-1)*Bcolwidth+1 : i*Bcolwidth) = Lshift * B(: , (i-2)*Bcolwidth+1 : (i-1)*Bcolwidth);
end

% C: matrix that projects lifted state into reference trajectory space
% based on the relationship: z_norm(k) = obj.model.C * h(k);
%                            \tilde{y(k)} = paramMapSSM * scale_up * z_norm(k);
%                            y(k) = refIso * \tilde{y(k)};
% the matrix that projects lifted state into reference trajectory space in
% the SSMK case is y(k) = refIso * paramMapSSM * obj.model.C * h(k);
% C_SSMK = refIso * paramMapSSM * model.C;
C = kron( speye(MPC.horizon+1) , model.C);
nproj = size( model.C , 1 );

% Q: Error magnitude penalty
Q = kron( speye(MPC.horizon+1) , var.costRun); % error magnitude penalty (running cost) (default 0.1)
Q(end-nproj+1 : end , end-nproj+1 : end) = var.costTer;    % (terminal cost) (default 100)
% Q = kron( speye(obj.horizon+1) , diag([obj.cost_running, 0, obj.cost_running]) ); % error magnitude penalty (running cost) (default 0.1)
% Q(end-nproj+1 : end , end-nproj+1 : end) = diag([obj.cost_terminal, 0, obj.cost_terminal]);    % (terminal cost) (default 100)

% R: Input magnitude penalty
R = kron( speye(MPC.horizon) , eye(model.params.m) * var.costInp );  % input magnitude penalty (for flaccy use 0.5e-2) (new videos used 0.5e-3)

% H, G, D
H = B' * C' * Q * C * B + R;
G = 2 * A' * C' * Q * C * B;
D = -2 * Q * C * B;

% set outputs
cost.H = H; cost.G = G; cost.D = D; % constructed matrices
cost.A = A; cost.B = B; cost.C = C; cost.Q = Q; cost.R = R; % component matrices
H = full(H);
G = full(G);
D = full(D);

% num = 2*model.params.m;     % number of input bound constraints

% F: input_bounds

F = []; E = [];     % initialize empty matrices
c = [];

if ~isempty( var.Inb ) && size( var.Inb , 1 ) ~= model.params.m
    var.Inb = kron( ones( model.params.m , 1 ) , var.Inb );
end

if ~isempty( var.Inb )
    num = 2*model.params.m;       % number of input bound constraints
    
    % F: input_bounds
    Fbounds_i = [ -speye(model.params.m) ; speye(model.params.m) ];    % diagonal element of F, for bounded inputs
    Fbounds = sparse( num * (MPC.horizon+1) , size(cost.B,2) );  % no constraints, all zeros
    Fbounds( 1:num*MPC.horizon , 1:MPC.horizon*model.params.m ) = kron( speye(MPC.horizon) , Fbounds_i );     % fill in nonzeros
    F = [ F ; Fbounds ];    % append matrix
    
    % E: input_bounds (just zeros)
    Ebounds = sparse( num * (MPC.horizon+1) , size(cost.B,1) );  % no constraints, all zeros
    E = [ E ; Ebounds ];    % append matrix
    
    % c: input_bounds
    input_bounds_sc = scaledown_u( var.Inb' )';   % scaled down the input bounds

    cbounds_i = [ -input_bounds_sc(:,1) ; input_bounds_sc(:,2) ]; % [ -umin ; umax ]
    cbounds = zeros( num * (MPC.horizon+1) , 1);    % initialization
    cbounds(1 : num*MPC.horizon) = kron( ones( MPC.horizon , 1 ) , cbounds_i );     % fill in nonzeros
    c = [ c ; cbounds ];    % append vector
end

% input_slopeConst
if ~isempty( var.InSlp)
    % F: input_slopeConst
    Fslope_i = speye(model.params.m);
    Fslope_neg = [ kron( speye(MPC.horizon-1) , -Fslope_i ) , sparse( model.params.m * (MPC.horizon-1) , model.params.m ) ];
    Fslope_pos = [ sparse( model.params.m * (MPC.horizon-1) , model.params.m ) , kron( speye(MPC.horizon-1) , Fslope_i ) ];
    Fslope_top = Fslope_neg + Fslope_pos;
    Fslope = [ Fslope_top ; -Fslope_top];
    F = [ F ; Fslope ];     % append matrix
    
    % E: input_slopeConst (just zeros)
    E = [ E ; sparse( 2 * model.params.m * (MPC.horizon-1) , size(cost.B,1) ) ];
    
    % c: input_slopeConst
    
    slope_lim = var.InSlp * mean( model.params.scale.u_factor );  % scale down the 2nd deriv. limit
    
    cslope_top = slope_lim * ones( model.params.m * (MPC.horizon-1) , 1 );
    cslope = [ cslope_top ; cslope_top ];
    c = [ c ; cslope ];     % append vector
end

% input_smoothConst
if ~isempty( var.InSmt )
    % F: input_smoothConst
    Fsmooth_i = speye(model.params.m);
    Fsmooth_lI = [ kron( speye(MPC.horizon-2) , Fsmooth_i ) , sparse( model.params.m * (MPC.horizon-2) , 2 * model.params.m ) ];
    Fsmooth_2I = [ sparse( model.params.m * (MPC.horizon-2) , model.params.m ) , kron( speye(MPC.horizon-2) , -2*Fslope_i ) , sparse( model.params.m * (MPC.horizon-2) , model.params.m ) ];
    Fsmooth_rI = [ sparse( model.params.m * (MPC.horizon-2) , 2 * model.params.m ) , kron( speye(MPC.horizon-2) , Fslope_i ) ];
    Fsmooth_top = Fsmooth_lI + Fsmooth_2I + Fsmooth_rI;
    Fsmooth = [ Fsmooth_top ; -Fsmooth_top ];
    F = [ F ; Fsmooth ];
    
    % E: input_smoothConst
    E = [ E ; sparse( 2 * model.params.m * (MPC.horizon-2) , size(cost.B,1) ) ];
    
    % c: input_smoothConst
    smooth_lim = model.params.Ts^2 * var.InSmt * mean( model.params.scale.u_factor );  % scale down the 2nd deriv. limit
    csmooth = smooth_lim * ones( size(Fsmooth,1) ,1);
    c = [ c ; csmooth ];
end


% set outputs

c = c;
L = F + E * cost.B;
M = E * cost.A;
L = full(L);
M = full(M);
return
end