function [Param, Flag, Opt] = plasmaAlphaExp()

    % Simulation settings
    Flag.model = "plasma";               % Model selection
    Flag.verbose = 4;
   
    Flag.scheme = "coupled";            % "coupled"/"split"                                                 (scheme.m)
    Flag.method = "fsolve";             % "fsolve"/"newton"                                                 (scheme.m)
    
    Flag.adaptive = true;               % fixed time step or adaptive                                       (simulate.m)
    Flag.CheckGradients = false;        % checks for jacobian correctness                                   (scheme.m)
        
    Flag.genterm = 'non-const';         % non-const/const;               const is mainly for testing        (scheme.m)
    Flag.compAlpha = true;              % whether or not to compute alpha (or beta) at the end              (postProcess.m)
    Flag.alpha = "exp";                 % "const"/"exp"                  controls how the alpha is treated  (scheme.m)
    
   
    Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)
    
    Flag.loadSol = "alphaExpBeta7e5-100sec";        % "no"/"checkpoint"/"SaveName"                                                (simulate.m)
    Flag.saveSol = "alphaExpBeta50kVn2";              % "no"/"SaveName"                                                   (postProcess.m)


    % Hyperparameters 
    Param.lr = 51;                   % Space grid points
    Param.K = 1000;                  % Time grid points
    Param.dt = 1e-10;                 % Time separation [s]
    Param.T = 100;                   % Total simulation time [s]
    Param.V0 = 3.017466567237710e+04;                    % Voltage at r=1 and t=1 [V]
    Param.ZhengIdx = 54;
    
    Param.mup = 1.86e-4;             % Mobility for positive species
    Param.mun = 5e-2;                % Mobility for negative species
    
    Param.stepBuffer = 0.75;         % Buffer for adaptive time stepping
    Param.scalingPower = 0.25;       % Scaling power for dt adjustment
    Param.dtMin = 1e-15;             % Minimum allowed timestep
    Param.dtMax = 1e-2;              % Maximum allowed timestep
       
    % Fsolve options
    Opt.Display = "off";                         % Display option: "off", "iter", "final", etc.
    Opt.SpecifyObjectiveGradient = false;        % Use automatic finite difference for gradient?
end
