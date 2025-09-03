function [B,Br,regErrorsNorm,regErrorsAvg] = fitControlMatricesCV(yuxData,RDInfoCtrl,IMInfoCtrl, regBr)
% Assume graph type approach where reduced coordinates are a linear
% projection. xuTraj is a list of full trajectories (of potentially delay 
% embedded data) where first column is time, second are the embedded
% observables, third the control and last the reduced coordinates 

redynFun = RDInfoCtrl.reducedDynamics.map;
monomials = RDInfoCtrl.reducedDynamics.phi;

t = [];       % time values
Xr = [];      % reduced coordinates at time k
U = [];       % controls at time k
dDvCdt = [];  % time derivatives at time k, difference between  
              % observables and autonomous parametrization 
dXrdt = [];   % time derivatives of reduced coordinates at time k
apprxOrd = 2; % approximation order for the derivative
% Data in matrices
for ii = 1:size(yuxData,1)
    t_in = yuxData{ii,1};
    Xr_in = yuxData{ii,2}; 
    U_in = yuxData{ii,3}; 

    [dXridt,Xri,ti] = finiteTimeDifferenceCV(Xr_in,t_in,apprxOrd);
    Ui = U_in(:, 1+apprxOrd:end-apprxOrd);
    t = [t ti]; Xr = [Xr Xri]; dXrdt = [dXrdt dXridt]; 
    U = [U Ui]; 
end
reDyn = redynFun(Xr);
reDynMono = monomials(Xr);
% Fit reduced dynamics control matrix
% reshape 
deltaDerivatives = dXrdt - reDyn;
% Feature construction (expand the Kronecker product of each column by column)
n_samples = size(U, 2);
ZkronU = zeros(size(reDynMono,1)*size(U,1), n_samples);

for i = 1:n_samples
    ZkronU(:, i) = kron(reDynMono(:, i), U(:, i));
end

X = [U; ZkronU];  % (n_samples x total_features)
Y = deltaDerivatives;            % (n_samples x output_dim)
% Learn whole B matrix
[Bapd,~,~] = ridgeRegression(X, Y, 0.01*ones(size(t)), [], regBr);
B  = Bapd(:, 1: size(U,1));                   % size: n_out x n_u
Br = Bapd(:, size(U,1)+ 1 : end); % reshape if needed

regErrorNorm = mean(sqrt(sum((deltaDerivatives - B*U - Br*ZkronU).^2)))/...
             max(sqrt(sum((deltaDerivatives).^2)));

regErrorAvg = mean(sqrt(sum((deltaDerivatives - B*U - Br*ZkronU).^2)));
               
% Output errors of the regression in %
regErrorsNorm = regErrorNorm*100;
regErrorsAvg = regErrorAvg*100;
end