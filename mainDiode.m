clear all;

addpath Utils
addpath init
addpath src
addpath src\assemblers\

% Hyperparameters Diode
Param = struct();
Param.K = 1;                   %Time grid points
Param.lr = 101;
Param.dt = 1e-14;               % Time separation  [s]
% Param.T  = 1e-6;
% Param.V0 = 0;                   % Starting voltage (r=1) [V]
Param.VT = 1.4;                 % Ending voltage (r=end) [V]
Param.case = 3;
Param.tolError = 1e-3;


Param.stepBuffer = 0.8;
Param.scalingPower = 0.3;
Param.dtMin = 1e-17;
Param.dtMax = 1e-2;


% simulation settings
Flag.model = "diode";
Flag.scheme = "coupled";

Flag.verbose = true;
Flag.Nverbose = false;

Flag.method = "fsolve";         % "fsolve"/"newton"
Flag.adaptive = true;        
Flag.CheckGradients = false;

Flag.mesh = "linear";           % "linear"/"tanh"
Flag.VT = "linear";          % "linear"/"piecewise"

Flag.loadSol = "no";            % "no"/"SaveName"
Flag.saveSol = "DiodeCase3";            % "no"/"SaveName"


% Fsolve flags
Opt = optimoptions('fsolve'); 
Opt.Display = "none";                   % "off"/"iter"/"final"/"final-detailed"     
Opt.SpecifyObjectiveGradient = false;    % Jacobian or no Jacobian
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
Results.ASol = solve(ADati, Flag, Opt);

% %% Postprocessing
% Results = postProcess(Dati, ADati, Results, Flag);
% 
% %% Plot
% % Plotting Flags
% Flag.concentrationPlot = "all";    % "all"/"last"/"none"
% Flag.potentialPlot = "none";
% Flag.currentPlot = "all";
% 
% plotter(Results,Dati,Flag);
