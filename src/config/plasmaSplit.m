function [Param, Flag, Opt] = plasmaSplit()

    % Simulation settings
    Flag.model = "plasma";               % Model selection
    Flag.verbose = 2;
   
    Flag.scheme = "split";            % "coupled"/"split"                                                 (scheme.m)
    Flag.method = "fsolve";             % "fsolve"/"newton"                                                 (scheme.m)
    
    Flag.adaptive = true;               % fixed time step or adaptive                                       (simulate.m)
    Flag.CheckGradients = false;        % checks for jacobian correctness                                   (scheme.m)
        
    Flag.genterm = 'non-const';         % non-const/const;               const is mainly for testing        (scheme.m)
    Flag.compAlpha = true;              % whether or not to compute alpha (or beta) at the end              (postProcess.m)
    Flag.alpha = "exp";                 % "const"/"exp"                  controls how the alpha is treated  (scheme.m)
    Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)
    
    Flag.loadSol = "no";        % "no"/"checkpoint"/"SaveName"                                                (simulate.m)
    Flag.saveSol = "temp";              % "no"/"SaveName"      


    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "last";
    
    % Hyperparameters 
    Param.lr = 51;                   % Space grid points
    Param.K = 1000;                  % Time grid points
    Param.dt = 1e-4;                 % Time separation [s]
    Param.T = 100;                   % Total simulation time [s]
    Param.V0 = 0;                    % Voltage at r=1 and t=1 [V]
    Param.tolError = 1e-3;          % Adaptive time stepping tollerance

    Param.mup = 1.86e-4;             % Mobility for positive species
    Param.mun = 5e-2;                % Mobility for negative species
    
    Param.stepBuffer = 0.8;         % Buffer for adaptive time stepping
    Param.scalingPower = 0.3;       % Scaling power for dt adjustment
    Param.dtMin = 1e-15;             % Minimum allowed timestep
    Param.dtMax = 5e-4;              % Maximum allowed timestep
       
    % Fsolve options
    Opt.Display = "off";                         % Display option: "off", "iter", "final", etc.
    Opt.SpecifyObjectiveGradient = false;        % Use automatic finite difference for gradient?

end



