clear all;

addpath utils
addpath init 
addpath src

% Hyper parameters Plasma
Param = struct();
Param.lr = 101;
Param.K = 1000;                 % Time grid points
Param.dt = 1e-4;                % Time separation  [s]
Param.T = 2;                    % Total simulation time
Param.V0 = 0;                   % Voltage at r=1 and t=1  [V]
% Param.VT = 0;                 % Ending voltage at r=1 and t=K*dt  [V]
Param.tolError = 1e-3;          % Adaptive time stepping tollerance
% Param.alpha = 11396.798;      % for alpha constant    
% Param.beta = 342439.622;      % for alpha exponential

% Controls the time scaling: scale = buffer*(tolError / relativeError)^scalePow;
Param.stepBuffer = 0.75;  % lower to reduce swinging in dt
Param.scalingPower = 0.25; % lower to reduce swinging in dt  

% Bounds dt
Param.dtMin = 1e-15;
Param.dtMax = 1e-2;



% Simulation settings
Flag.model = "plasma";


Flag.method = "fsolve";             % "fsolve"/"newton"
Flag.verbose = true;
Flag.Nverbose = 0;

Flag.adaptive = true;        
Flag.CheckGradients = false;

Flag.genterm = 'non-const';         % non-const/const;
Flag.compAlpha = true;              % whether to compute alpha or beta at the end
Flag.alpha = "exp";                 % "const"/"exp"
Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"

Flag.loadSol = "no";                % "no"/"SaveName"
Flag.saveSol = "temp";              % "no"/"SaveName"

% Newton solver Parameters

% this whole shenaningans and with set_options can probably be replaced by
% Opt = optimoptions('fsolve'); and add the ones needed and pass it already
% created.

% Fsolve options
Opt.Display = "off";                        % "off"/"iter"/"final"/"final-detailed"
Opt.SpecifyObjectiveGradient = false; 
Opt.FiniteDifferenceType = "forward";       % "central"/"forward"
% Opt.Algorithm = "levenberg-marquardt";    % "levenberg-marquardt

% Opt.MaxIterations = 400;
% Opt.MaxFunctionEvaluations = 50e4;
% Opt.OptimalityTolerance = 1e-6;
Opt.StepTolerance = 1e-6;
Opt.FunctionTolerance = 1e-6;


% Takes care of Parameters, Intial condition, and adimensionalization
[Param, Dati, ADati] = initPlasma(Param, Flag);

%% Solve
% Start logging
logFile = 'run_log.txt';
if exist(logFile, 'file')
    diary off;
    delete(logFile); % Optional: delete previous log
end
diary(logFile);     % Start recording to file
diary on;

tic;
try
    Results.ASol = solve(ADati, Flag, Opt);
catch ME
    disp(ME)
    fprintf('ERROR: %s\n', ME.message);
    diary off;
    return
end

Results.elapsedTime = toc;

%% Postprocessing
Results = postProcess(Dati, ADati, Results, Flag);

%% Plot
% Plotting Flags
Flag.concentrationPlot = "last";    % "all"/"last"/"none"
Flag.potentialPlot = "none";
Flag.currentPlot = "all";

Flag.generationPlot = "none";

plotter(Results,Dati,Flag);



diary off;