clc
clear all
load SSMP_model_0814.mat;
size_x = 2;                   % number of states
size_u = py_data.params.input_dim;                   % number of control input
size_f = 2;                   % number of reference
dis = 110;                     % location of the work space, unit: mm
[size_C, size_Cbar]=generate_functions(size_x,size_u,size_f,py_data);


function [size_C, size_Cbar]=generate_functions(size_x,size_u,size_f,py_data)
%%-------------------------------------------------------------------------        
% Get users' functions for generating sub-functions of users' optimization 
% problems from the trained SSM-model in .\SSM Model Learning folder.
%
% Last update: April/07/2024 by Yongjun Yan
% Revision history: The reference is changed as the trajectory at the SSM
%  x     : The normalized state at the low-dimensional coordinates.
%  u     : The control variables.
%  z     : The observables without delay-embedding.
%  r     : The reference.
syms x u r z
for n = 1:1:size_x
    syms (sprintf('x%01d',n));
    x(n) = sym(sprintf('x%01d',n));
end
for m = 1:1:size_u
    syms (sprintf('u%01d',m));
    u(m) = sym(sprintf('u%01d',m));
end
for l = 1:1:size_f
    syms (sprintf('r%01d',l));
    r(l) = sym(sprintf('r%01d',l));
end
for o = 1:1:7
    syms (sprintf('z%01d',o));
    z(o) = sym(sprintf('z%01d',o));
end

%% ------------------------------------------------------------------------
% Do not modify below this line 
%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------        
% Get users' functions for generating sub-functions of users' optimization
%  problems to evaluate values of function and costraints values 
%  and to compute derivatives.
%
% Last update: Feb/21/2023 by Yongjun Yan
%-------------------------------------------------------------------------
% Users should provide following 5 functions in the symbolic form.
% Phi:    terminal cost function
% L:      incremental cost function
% f:      equations of system dynamics 
% C:      input and state mixed constraints
% Cbar:   state only constraints
%--------------------------------------------------------------------------
% Example:
% if your cost function 
%  J=0.5*y(N)'*S*y(N) + 0.5*\Sum_(k=0)^(k=N-1)(y(k)'*Q*y(k) + y(k)'*R*y(k),
% and your eqautaion of system dynamics x(1,k+1) = x(1,k)+Ts*u(1,k), 
%                                       x(2,k+1) = x(2,k)+Ts*x(1,k)*u(2,k),
%                       and constraints x(1,k) + u(1,k) <= -1,       
%                                       x(2,k) + u(2,k) >= 5,
%                                       -2 <= x(1,k) <= 2, 
%                                       x(1,k)^2 + x(2,k)^2 <= 4.
%
% Then, first define all parameters and define the functions based on
% the rules of Symbolic Math Toolbox
% Phi, L, f, C, and Cbar will be
%
% % define parameters
% s = 1e4;  q = 1,  r= 1e3; 
%
% % Terminal cost
% Phi = 0.5*s*(x1^2+x2^2);
%
% % Incremental cost
% L = 0.5*q*(x1^2+x2^2) + 0.5*r*(u1^2+u2^2);
%
% % Equations of system dynamics
% f = [x1+Ts*u1;
%      x2+Ts*x1*u2];
% 
% % Input and state mixed constriants
% C = [ x1+u1+1;
%       -x2-u2+5];
%        
% % State only constraints 
% Cbar = [ -x1-2;
%           x1-2;
%           x1^2+x2^2-4];
%--------------------------------------------------------------------------     
%% ------------------------------------------------------------------------
% Users' functions: Please provide your functions
% -------------------------------------------------------------------------
%
% The x is the state vector composed with the Cartesian position of the two
% feature points, after minusing the equilibrium point;
% 6Dof = [6Dof_x, 6Dof_y, 6Dof_z] is the Cartesian position get from the
% 6Dof sensor;
% Sec = [Sec_x, Sec_y, Sec_z] is the Cartesian position derived by
% projecting the Cartesian position of the 6Dof sensor along the z-axis
% direction, distance is 1 mm away from the 6Dof sensor;
% x = [Sec, 6Dof]';



% Equations of system dynamics
Dx_lowDim = py_data.model.R(x',u');
% f_Org = diag([0.8,1,0.4])*x' + py_data.model.Ts*Dx_lowDim_Sim;
f = x' + py_data.model.Ts*Dx_lowDim;
% f = vpa(f_Org,3); 
x_su = py_data.scale.scaleup.y(x);
y = py_data.model.IMInfoCtrl.parametrization.map(x_su');



% define the model parameters

pf = diag([100,100]);

% Penalty Costs
q = diag([1,1]);

% Terminal cost
Phi = 0.5*(  (y(2)-r1)*pf(1,1)*(y(2)-r1) + (y(3)-r2)*pf(2,2)*(y(3)-r2) );
% Phi = vpa(Phi_Org,3);

% Incremental cost

% L_Org = 0.5*(   (y1-ref1)*q(1,1)*(y1-ref1) + (y1-ref1)*q(1,2)*(y2-ref2) ...
%               + (y2-ref2)*q(2,1)*(y1-ref1) + (y2-ref2)*q(2,2)*(y2-ref2)     );
L = 0.5*(  (y(2)-r1)*q(1,1)*(y(2)-r1) + (y(3)-r2)*q(2,2)*(y(3)-r2) );
% L = vpa(L_Org,3);

u_ori = py_data.scale.scaleup.u(u);
% Input and state mixed constriants
C =  [u_ori(1)-1100;
      u_ori(2)-1100;
     -u_ori(1)-1100;
     -u_ori(2)-1100];

% State only constraints 
Cbar = [x1-1;
        x2-1;
       -x1-1;
       -x2-1];
%%
scaleup_u = py_data.scale.scaleup.u(u);
scaleup_y = py_data.scale.scaleup.y(x);
scaledown_u = py_data.scale.scaledown.u(u);
scaledown_y = py_data.scale.scaledown.y(x);
% chart map
chart = py_data.model.IMInfoCtrl.chart.map(z'); % Project down
% parameterization map
param = py_data.model.IMInfoCtrl.parametrization.map(x');
%% ------------------------------------------------------------------------
% Do not modify below this line 
%-------------------------------------------------------------------------

if isempty(Phi)
    Phi = zeros(0,1);
    matlabFunction(Phi, 'file', 'valPhi','vars',[x, u, r]);
else
    matlabFunction(Phi, 'file', 'valPhi','vars',{x, u, r}); 
end

if isempty(L)
    L = zeros(0,1);
    matlabFunction(L, 'file', 'valL','vars',[x, u, r]); 
else
    matlabFunction(L, 'file', 'valL','vars',{x, u, r}); 
end

if isempty(Dx_lowDim)
    Phi = zeros(0,1);
    matlabFunction(Dx_lowDim, 'file', 'valDx','vars',[x, u]);
else
    matlabFunction(Dx_lowDim, 'file', 'valDx','vars',{x, u}); 
end

if isempty(f)
    f = zeros(0,1);
    matlabFunction(f, 'file', 'valf','vars',[x, u]); 
else
    matlabFunction(f, 'file', 'valf','vars',{x, u}); 
end

if isempty(C)
    C = zeros(0,1);
    matlabFunction(C, 'file', 'valC','vars',[x, u]);   
else
    matlabFunction(C, 'file', 'valC','vars',{x, u});
end

if isempty(Cbar)
    Cbar = zeros(0,1);
    matlabFunction(Cbar, 'file', 'valCbar','vars',x); 
else
    matlabFunction(Cbar, 'file', 'valCbar','vars',{x}); 
end

if isempty(scaleup_u)
    scaleup_u = zeros(0,1);
    matlabFunction(scaleup_u, 'file', 'scaleup_u','vars',[u]);
else
    matlabFunction(scaleup_u, 'file', 'scaleup_u','vars',{u}); 
end

if isempty(scaleup_y)
    scaleup_y = zeros(0,1);
    matlabFunction(scaleup_y, 'file', 'scaleup_y','vars',[x]);
else
    matlabFunction(scaleup_y, 'file', 'scaleup_y','vars',{x}); 
end

if isempty(scaledown_u)
    scaledown_u = zeros(0,1);
    matlabFunction(scaledown_u, 'file', 'scaledown_u','vars',[u]);
else
    matlabFunction(scaledown_u, 'file', 'scaledown_u','vars',{u}); 
end

if isempty(scaledown_y)
    scaledown_y = zeros(0,1);
    matlabFunction(scaledown_y, 'file', 'scaledown_y','vars',[x]);
else
    matlabFunction(scaledown_y, 'file', 'scaledown_y','vars',{x}); 
end

if isempty(chart)
    chart = zeros(0,1);
    matlabFunction(chart, 'file', 'chart','vars',[z]);
else
    matlabFunction(chart, 'file', 'chart','vars',{z}); 
end

if isempty(param)
    param = zeros(0,1);
    matlabFunction(param, 'file', 'param','vars',[x]);
else
    matlabFunction(param, 'file', 'param','vars',{x}); 
end

%%
if isempty(Phi)
    Phix = zeros(0,1);
    Phixx = zeros(0,1);
    
    matlabFunction(Phix, 'file', 'DPhiDx','vars',x);
    matlabFunction(Phixx, 'file', 'D2PhiDxx','vars',x); 
else
    for j=1:size_x
        Phix(1,j) = diff(Phi,(sprintf('x%1d',j)));
    end
    matlabFunction(Phix, 'file', 'DPhiDx','vars',{x, u, r});   

    for j = 1:size_x
        for k = 1:size_x
            Phixx(j,k) = diff(Phix(1,k),(sprintf('x%1d',j)));
        end
    end
    matlabFunction(Phixx, 'file', 'D2PhiDxx','vars',{x, u, r}); 
end

%%

if isempty(L)
    Lx = zeros(0,1);
    Lu = zeros(0,1);
    Lxx = zeros(0,1);
    Lxu = zeros(0,1);
    Luu = zeros(0,1);
    
    matlabFunction(Lx, 'file', 'DLDx','vars',[x, u]);
    matlabFunction(Lu, 'file', 'DLDu','vars',[x, u]);
    matlabFunction(Lxx, 'file', 'D2LDxx','vars',[x, u]);
    matlabFunction(Lxu, 'file', 'D2LDxu','vars',[x, u]); 
    matlabFunction(Luu, 'file', 'D2LDuu','vars',[x, u]);
else
    for j=1:size_x
        Lx(1,j) = diff(L,(sprintf('x%1d',j)));
    end
    matlabFunction(Lx, 'file', 'DLDx','vars',{x, u, r});   

    for j=1:size_u
        Lu(1,j) = diff(L,(sprintf('u%1d',j)));
    end
    matlabFunction(Lu, 'file', 'DLDu','vars',{x, u, r});   

    for j = 1:size_x
        for k = 1:size_x
            Lxx(j,k) = diff(Lx(1,k),(sprintf('x%1d',j)));
        end
    end
    matlabFunction(Lxx, 'file', 'D2LDxx','vars',{x, u, r}); 

    for j = 1:size_u
        for k = 1:size_x
            Lxu(j,k) = diff(Lx(1,k),(sprintf('u%1d',j)));
        end
    end
    matlabFunction(Lxu, 'file', 'D2LDxu','vars',{x, u, r}); 

    for j = 1:size_u
        for k = 1:size_u
            Luu(j,k) = diff(Lu(1,k),(sprintf('u%1d',j)));
        end
    end
    matlabFunction(Luu, 'file', 'D2LDuu','vars',{x, u, r}); 

end


%%
size_f = size(f,1);
save size_f.mat size_f;
if isempty(f)
    fx = zeros(0,1);
    fu = zeros(0,1);
    fxx = zeros(0,1);
    fxu = zeros(0,1);
    fuu = zeros(0,1);
    
    matlabFunction(fx, 'file', 'DfDx','vars',[x, u]);
    matlabFunction(fu, 'file', 'DfDu','vars',[x, u]);
    matlabFunction(fxx, 'file', 'D2fDxx','vars',[x, u]);
    matlabFunction(fxu, 'file', 'D2fDxu','vars',[x, u]); 
    matlabFunction(fuu, 'file', 'D2fDuu','vars',[x, u]);
else
    for j=1:size_x
        for k = 1:size_x
            fx(j,k) = diff(f(j),(sprintf('x%1d',k)));
        end 
    end
    matlabFunction(fx, 'file', 'DfDx','vars',{x, u});   

    for j=1:size_x
        for k = 1:size_u
            fu(j,k) = diff(f(j),(sprintf('u%1d',k)));
        end 
    end
    matlabFunction(fu, 'file', 'DfDu','vars',{x, u});   

    for n=1:size_x
        for j = 1:size_x
            for k = 1:size_x
                fxx(j,k,n) = diff(fx(n,k),(sprintf('x%1d',j)));
            end
        end
    end
    matlabFunction(fxx, 'file', 'D2fDxx','vars',{x, u}); 

    for n=1:size_x
        for j = 1:size_u
            for k = 1:size_x
                fxu(j,k,n) = diff(fx(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(fxu, 'file', 'D2fDxu','vars',{x, u}); 

    for n=1:size_x
        for j = 1:size_u
            for k = 1:size_u
                fuu(j,k,n) = diff(fu(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(fuu, 'file', 'D2fDuu','vars',{x, u}); 
end


%%
size_C = size(C,1);
save size_C.mat size_C;
if size_C == 0
    Cx = zeros(0,1);
    Cu = zeros(0,1);
    Cxx = zeros(0,1);
    Cxu = zeros(0,1);
    Cuu = zeros(0,1);
    
    matlabFunction(Cx, 'file', 'DCDx','vars',[x, u]);
    matlabFunction(Cu, 'file', 'DCDu','vars',[x, u]);
    matlabFunction(Cxx, 'file', 'D2CDxx','vars',[x, u]);
    matlabFunction(Cxu, 'file', 'D2CDxu','vars',[x, u]); 
    matlabFunction(Cuu, 'file', 'D2CDuu','vars',[x, u]);
else
    for j=1:size_C
        for k = 1:size_x
            Cx(j,k) = diff(C(j),(sprintf('x%1d',k)));
        end
    end
    matlabFunction(Cx, 'file', 'DCDx','vars',{x, u});   

    for j=1:size_C
        for k = 1:size_u
            Cu(j,k) = diff(C(j),(sprintf('u%1d',k)));
        end
    end
    matlabFunction(Cu, 'file', 'DCDu','vars',{x, u});   

    for n=1:size_C
        for j = 1:size_x
            for k = 1:size_x
                Cxx(j,k,n) = diff(Cx(n,k),(sprintf('x%1d',j)));
            end
        end
    end
    matlabFunction(Cxx, 'file', 'D2CDxx','vars',{x, u}); 

    for n=1:size_C
        for j = 1:size_u
            for k = 1:size_x
                Cxu(j,k,n) = diff(Cx(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(Cxu, 'file', 'D2CDxu','vars',{x, u}); 

    for n=1:size_C
        for j = 1:size_u
            for k = 1:size_u
                Cuu(j,k,n) = diff(Cu(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(Cuu, 'file', 'D2CDuu','vars',{x, u}); 
end

%%
size_Cbar = size(Cbar,1);
save size_Cbar.mat size_Cbar;
if size_Cbar == 0
   Cbarx = zeros(0,1);
   Cbarxx = zeros(0,1);
   matlabFunction(Cbarx, 'file', 'DCbarDx','vars',x);
   matlabFunction(Cbarxx, 'file', 'D2CbarDxx','vars',x); 
else
    for j=1:size_Cbar
        for k = 1:size_x
            Cbarx(j,k) = diff(Cbar(j),(sprintf('x%1d',k)));
        end
    end
    matlabFunction(Cbarx, 'file', 'DCbarDx','vars',{x});   

    for n=1:size_Cbar
        for j = 1:size_x
            for k = 1:size_x
                Cbarxx(j,k,n) = diff(Cbarx(n,k),(sprintf('x%1d',j)));
            end
        end
    end
    matlabFunction(Cbarxx, 'file', 'D2CbarDxx','vars',{x}); 
end

return
end