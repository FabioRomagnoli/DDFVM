1function [Param, Flag, Opt] = plasmaAlphaConst()
    % JACOBIAN DOESN'T WORK FOR THIS YET
    % This does start running but suuper slow
    
    % simulation settings
    Flag.model = "plasma";               % Model selection
    Flag.scheme = "coupled";
    Flag.verbose = true;

    Flag.scheme = "coupled";            % "coupled"/"split"                                                 (scheme.m)
    Flag.method = "fsolve";             % "fsolve"/"newton"                                                 (scheme.m)
    
    Flag.adaptive = true;               % fixed time step or adaptive                                       (simulate.m)
    Flag.CheckGradients = false;        % checks for jacobian correctness                                   (scheme.m)
        
    Flag.genterm = 'non-const';         % non-const/const;
    Flag.alpha = "const";                 % "const"/"exp"
    Flag.compAlpha = true;              % whether or not to compute alpha (or beta) at the end              (postProcess.m)


    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "last";

    Flag.loadSol = "no";                     % "no"/"checkpoint"/"SaveName"                                                (simulate.m)
    Flag.saveSol = "plasmaConstAlpha";       % "no"/"SaveName"                                                   (postProcess.m)


    % Hyperparameters
    Param.lr = 101;                   % Space grid points
    Param.K = 1000;                  % Time grid points
    Param.dt = 1e-4;                 % Time separation [s]
     
    Param.T = 1;                    % Total simulation time
    Param.V0 = 0;                   % Voltage at r=1 and t=1  [V]
    % Param.alpha = 11396.798;      % for alpha constant 

    % Maybe different alphas? This does run, maybe go to the  end?
    Param.alpha = 100;

    % Fsolve flags
    Opt.Display = "off";                         % Display option: "off", "iter", "final", etc.
    Opt.SpecifyObjectiveGradient = false;        % Use automatic finite difference for gradient?
end