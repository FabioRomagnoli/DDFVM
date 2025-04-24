function [Param, Flag, Opt] = diodeCase3()

    % simulation settings
    Flag.model = "diode";
    Flag.scheme = "coupled";
    Flag.verbose = true;
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.adaptive = true;        
    Flag.mesh = "linear";           % "linear"/"tanh"

    Flag.loadSol =  ".\diode\diodeCase3";
    Flag.saveSol = ".\diode\diodeCase3";
        
    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "last";
    


    % Hyperparameters Diode
    Param = struct();
    Param.K = 50;                   %Time grid points
    Param.lr = 101;
    Param.dt = 1e-10;               % Time separation  [s]
    Param.VT = 1.4;                 % Ending voltage (r=end) [V]
    Param.case = 3;
    
    Param.stepBuffer = 0.8;
    Param.scalingPower = 0.3;
    Param.dtMin = 1e-17;
    Param.dtMax = 1e-2;
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian
end



