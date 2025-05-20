function [Param, Flag, Opt] = diodeCase3()

    % simulation settings

    Flag.model = "diode";
    Flag.scheme = "coupled";
    Flag.verbose = 3;
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.adaptive = false;        
    Flag.mesh = "linear";           % "linear"/"tanh"

    Flag.loadSol =  "no";
    Flag.saveSol = ".\diode\diodeConv51NotAdaptive";
        
    % plotting
    Flag.concentrationPlot = "last";    % Options: "all", "last", "none"
    Flag.potentialPlot = "last";
    Flag.currentPlot = "last";
    


    % Hyperparameters Diode
    Param = struct();
    Param.K = 50;                   %Time grid points

    Param.lr = 51;
    % Param.lr = 101;
    % Param.lr = 201;
    % Param.lr = 401;
    % Param.lr = 801;

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



