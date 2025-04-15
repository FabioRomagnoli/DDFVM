function [Param, Flag, Opt] = diodeCase1()

    % simulation settings
    Flag.model = "diode";
    Flag.scheme = "coupled";
    Flag.verbose = false;
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.adaptive = false;        
    Flag.mesh = "linear";           % "linear"/"tanh"
    Flag.saveSol = ".\diode\diodeCase1";
        
    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "none";
    
    % Hyperparameters Diode
    Param.K = 1;                   %Time grid points
    Param.lr = 101;
    Param.dt = 1e-10;               % Time separation  [s]
    Param.V0 = 1;
    Param.VT = 1;                 % Ending voltage (r=end) [V]
    Param.case = 1;
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian
end
