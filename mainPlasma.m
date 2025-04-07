clear all;

addpath utils
addpath plasma 

% Hyper parameters Plasma
Param = struct();
Param.lr = 101;
Param.K = 100;                  % Time grid points
Param.dt = 1e-4;                % Time separation  [s]
Param.T = 100 * Param.K*Param.dt;
Param.V0 = 0;                 % Voltage at r=1 and t=1  [V]
% Param.VT = 0;                 % Ending voltage at r=1 and t=K*dt  [V]
Param.N = 1e7;                  % density constant [m-3]
Param.tolError = 1e-3;
Param.alpha = 11396.798;  
Param.Beta = 342439.622;

Param.stepBuffer = 0.8;
Param.scalingPower = 0.3;
Param.dtMin = 1e-11;
Param.dtMax = 1e-2;


% Simulation settings
Flag.method = "fsolve";         % "fsolve"/"newton"
Flag.adaptive = true ;        
Flag.CheckGradients = false;

Flag.genterm = 'non-const';         % non-const/const;
Flag.compAlpha = true;
Flag.alpha = "exp";         % "const"/"exp"
Flag.VT = "linear";             % "linear"/"piecewise"

Flag.saveSol = "no";            % "no"/"SaveName"
Flag.loadSol = "no";            % "no"/"SaveName"


% Fsolve options
Opt.Display = "iter-detailed";                    % "off"/"iter"/"final"/"final-detailed"
Opt.SpecifyObjectiveGradient = false; 
Opt.FiniteDifferenceType = "forward";   % "central"/"forward"
% Opt.Algorithm = "trust-region-dogleg";
% Opt.MaxIterations = 400;
% Opt.MaxFunctionEvaluations = 40e4;
% Opt.OptimalityTolerance = 1e-6;
% Opt.StepTolerance = 1e-6;
% Opt.FunctionTolerance = 1e-6;


% Takes care of Parameters, Intial condition, and adimensionalization
[Param, Dati, ADati] = initPlasma(Param, Flag);

% Solve
ASol = solvePlasma(ADati, Flag, Opt);

%% Postprocessing
[Sol, Res] = postPlasma(ADati, Dati, ASol, Flag);

%% Plot
% Plotting Flags
Flag.concentrationPlot = "none";    % "all"/"last"/"none"
Flag.potentialPlot = "none";
Flag.currentPlot = "all";

Flag.generationPlot = "none";
Flag.compCurrentPlot = "none";

plotPlasma(Sol,Res,Dati,Flag);

