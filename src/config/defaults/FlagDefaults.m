function def = FlagDefaults()
    % SOLVER  
    def.scheme = "coupled";             % Coupling scheme: "coupled" or "split"
    def.method = "fsolve";              % Solver choice: "fsolve" or "newton"
    def.CheckGradients = false;         % Option to check the correctness of the Jacobian

    % VERBOSITY
    def.verbose = 2;                 % Verbosity flag for simulation output
    def.Nverbose = 0;                   % Verbosity level for the nonlinear solver

    % TIME 
    def.adaptive = true;                % Enable adaptive time stepping
    def.VT = "plateu";                  % Voltage evolution: "linear", "piecewise", or "plateu"
    def.EndVT = "plateu";                  % Voltage evolution: "linear", "piecewise", or "plateu"

    def.computeCurrents = true;               % Computes current at every time step

    % DIODE ONLY
    def.mesh = "tanh";           % "linear"/"tanh"

    % PLASMA
    def.genterm = 'non-const';          % Setting for the general term ('non-const' or 'const')
    def.compAlpha = true;               % Option to compute alpha (or beta) at postprocessing

    def.alpha = "exp";                  % Treatment of alpha: "const" or "exp"

    % SAVING AND LOADING
    def.loadSol = "no";         % Option for loading a previous solution: "no", "checkpoint", or a filename
    def.saveSol = "no";    


    % Plotting
    def.concentrationPlot = "none";    % Options: "all", "last", "none"
    def.potentialPlot = "none";
    def.currentPlot = "none";
    def.generationPlot = "none";
end