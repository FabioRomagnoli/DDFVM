function [Param, Flag, Opt] = diodeCase3()

    % simulation settings
    Flag.model = "diode";
    Flag.verbose = 2;

    Flag.scheme = "coupled";
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)

    Flag.adaptive = true;        
    Flag.mesh = "linear";           % "linear"/"tanh"

    Flag.loadSol =  "no";
    Flag.saveSol = ".\diode\diodeCase3inTime";
        
    

    % Hyperparameters Diode
    Param.K = 1000;                   %Time grid points
    Param.lr = 101;
    Param.T = 100;                   % Total simulation time [s]

    Param.dt = 1e-10;               % Time separation  [s]

    Param.V0 = 0;                    % Voltage at r=1 and t=1 [V]
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



