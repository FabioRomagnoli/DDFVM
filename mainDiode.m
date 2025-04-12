clear all;

addpath(genpath('.\init'))
addpath(genpath('.\utils'))
addpath(genpath('.\src'))



% Hyperparameters Diode
Param = struct();
Param.K = 100;                   %Time grid points
Param.lr = 101;
Param.dt = 1e-10;               % Time separation  [s]
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
Flag.Nverbose = 0;

Flag.method = "fsolve";         % "fsolve"/"newton"
Flag.adaptive = false;        
Flag.CheckGradients = true;

Flag.mesh = "linear";           % "linear"/"tanh"
Flag.VT = "linear";          % "linear"/"piecewise"

Flag.loadSol = "no";            % "no"/"SaveName"
Flag.saveSol = "DiodeCase3";            % "no"/"SaveName"


% Fsolve flags
Opt = optimoptions('fsolve'); 
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
Results.ASol = simulate(ADati, Flag, Opt);

%% Postprocessing
Results = postProcess(Dati, ADati, Results, Flag);

%% Plot
% Plotting Flags
Flag.concentrationPlot = "all";    % "all"/"last"/"none"
Flag.potentialPlot = "none";
Flag.currentPlot = "all";

plotter(Results,Dati,Flag);
