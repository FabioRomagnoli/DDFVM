clear all;

addpath Utils
addpath diode

% Hyperparameters Diode
Param = struct();
Param.K = 100;                   %Time grid points
Param.lr = 101;
Param.dt = 1e-10;               % Time separation  [s]
% Param.T  = 1e-6;
% Param.V0 = 0.7;                   % Starting voltage (r=1) [V]
Param.VT = 1.4;                 % Ending voltage (r=end) [V]
Param.case = 3;
Param.tol = 1e-3;
Param.Nverbose = 0;

% simulation settings
Flag.method = "newton";         % "fsolve"/"newton"
Flag.adaptive = true;        
Flag.CheckGradients = true;

Flag.mesh = "linear";           % "linear"/"tanh"
Flag.VT = "piecewise";          % "linear"/"piecewise"

Flag.saveSol = "no";            % "no"/"SaveName"
Flag.loadSol = "no";            % "no"/"SaveName"


% Fsolve flags
Opt.Display = "none";                   % "off"/"iter"/"final"/"final-detailed"     
Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian
Opt.Algorithm = "trust-region-dogleg";  % trust-region-dogleg
Opt.FiniteDifferenceType = "forward";   % "central"/"forward"
% Opt.MaxIterations = 400;
% Opt.MaxFunctionEvaluations = 29e3;
Opt.OptimalityTolerance = 1e-6;
Opt.StepTolerance = 1e-6;
Opt.FunctionTolerance = 1e-6;


% Takes care of Parameters, Intial condition, and adimensionalization
[Param, Dati, ADati] = initDiode(Param, Flag);

% Solve
ASol = solveDiode(ADati, Flag, Opt);

%% Postprocessing
[Sol, Res] = postDiode(ADati, Dati, ASol, Flag);

%% Plot
% Plotting Flags
Flag.concentrationPlot = "last";    % "all"/"last"/"none"
Flag.potentialPlot = "none";
Flag.currentPlot = "all";

plotDiode(Sol,Res,Dati,Flag);
