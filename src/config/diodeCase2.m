function [Param, Flag, Opt] = diodeCase2()

    % simulation settings
    Flag.model = "diode";
    Flag.scheme = "coupled";
    Flag.verbose = true;
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.adaptive = true;        
    Flag.mesh = "linear";           % "linear"/"tanh"
    Flag.saveSol = ".\diode\diodeCase2";
        
    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "none";
    

    % Hyperparameters Diode
    Param.K = 50;                   %Time grid points
    Param.lr = 101;
    Param.dt = 1e-10;               % Time separation  [s]
    Param.V0 = 1.4;
    Param.VT = 1.4;                 % Ending voltage (r=end) [V]
    Param.case = 2;
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian
end
