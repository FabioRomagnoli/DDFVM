function [Param, Flag, Opt] = plasmaConstGen()
    
    % simulation settings
    Flag.model = "plasma";               % Model selection
    Flag.verbose = true;
    Flag.Nverbose = 0;

    Flag.scheme = "coupled";            % "coupled"/"split"                                                 (scheme.m)
    Flag.method = "fsolve";             % "fsolve"/"newton"                                                 (scheme.m)
    
    Flag.adaptive = false;               % fixed time step or adaptive                                       (simulate.m)
    Flag.CheckGradients = false;        % checks for jacobian correctness                                   (scheme.m)
        
    Flag.genterm = 'const';             % non-const/const;               const is mainly for testing        (scheme.m)
    Flag.compAlpha = true;              % whether or not to compute alpha (or beta) at the end              (postProcess.m)
    Flag.alpha = "const";

    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "last";
  
    Flag.loadSol = "no";        % "no"/"checkpoint"/"SaveName"                                                (simulate.m)
    Flag.saveSol = "plasmaConstGen";              % "no"/"SaveName"                                                   (postProcess.m)

    % Hyperparameters
    Param.lr = 101;                   % Space grid points
    Param.K = 100;                  % Time grid points
    Param.dt = 1e-4;                 % Time separation [s]
     
    % Fsolve flags
    Opt.Display = "off";                         % Display option: "off", "iter", "final", etc.
    Opt.SpecifyObjectiveGradient = true;        % Use automatic finite difference for gradient?
end
