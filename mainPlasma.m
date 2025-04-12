clear all;

addpath(genpath('.\init'))
addpath(genpath('.\utils'))
addpath(genpath('.\src'))


% Hyper parameters Plasma
Param = struct();
Param.lr = 101;
Param.K = 1000;                 % Time grid points
Param.dt = 1e-3;                % Time separation  [s]
Param.T = 100;                    % Total simulation time
Param.V0 = 0;                   % Voltage at r=1 and t=1  [V]
% Param.VT = 0;                 % Ending voltage at r=1 and t=K*dt  [V]
Param.tolError = 1e-3;          % Adaptive time stepping tollerance
% Param.alpha = 11396.798;      % for alpha constant    
% Param.beta = 342439.622;      % for alpha exponential

% MUs
Param.mup = 1.86e-4;
Param.mun = 5e-2;


% Controls the time scaling: scale = buffer*(tolError / relativeError)^scalePow;
Param.stepBuffer = 0.75;  % lower to reduce swinging in dt
Param.scalingPower = 0.25; % lower to reduce swinging in dt  

% Bounds dt
Param.dtMin = 1e-15;
Param.dtMax = 1e-2;



% Simulation settings
Flag.model = "plasma";              
Flag.scheme = "coupled";            % "coupled"/"split"                                                 (scheme.m)

Flag.method = "fsolve";             % "fsolve"/"newton"                                                 (scheme.m)
Flag.verbose = true;                % verbosity of the output                                           (scheme.m)
Flag.Nverbose = 0;                  % verbsotiy of newton solver (comparable to opt.display)            (scheme.m)

Flag.adaptive = true;               % fixed time step or adaptive                                       (simulate.m)
Flag.CheckGradients = false;        % checks for jacobian correctness                                   (scheme.m)
    
Flag.genterm = 'non-const';         % non-const/const;               const is mainly for testing        (scheme.m)
Flag.compAlpha = true;              % whether or not to compute alpha (or beta) at the end              (postProcess.m)
Flag.alpha = "exp";                 % "const"/"exp"                  controls how the alpha is treated  (scheme.m)
Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)


Flag.loadSol = "checkpoint";        % "no"/"checkpoint"/"SaveName"                                                (simulate.m)
Flag.saveSol = "temp";              % "no"/"SaveName"                                                   (postProcess.m)



% Fsolve options
Opt = optimoptions('fsolve'); 
Opt.Display = "off";                        % "off"/"iter"/"final"/"final-detailed"
Opt.SpecifyObjectiveGradient = false;       % JACOBIAN
Opt.FiniteDifferenceType = "forward";       % "central"/"forward"
% Opt.MaxIterations = 400;
% Opt.MaxFunctionEvaluations = 4e4;
% Opt.OptimalityTolerance = 1e-6;
% Opt.StepTolerance = 1e-6;
% Opt.FunctionTolerance = 1e-6;

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
    % ACTUAL SOLVER
    Results.ASol = simulate(Dati, ADati, Flag, Opt);
catch ME
    fprintf(2, '\n%s\n', getReport(ME, 'extended'));
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