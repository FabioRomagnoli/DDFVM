function [Param, Flag, Opt] = diodeCase3()

    % simulation settings
    Flag.model = "diode";
    Flag.verbose = 1;

    Flag.scheme = "coupled";
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.VT = "linear";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)
    Flag.EndVT = "linear";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)

    Flag.adaptive = true;        
    Flag.mesh = "linear";           % "linear"/"tanh"

    Flag.loadSol =  "no";
    Flag.saveSol = ".\diode\diodeCase5Adaptive";
        
    
    % Hyperparameters Diode
    Param.K = 100;                   %Time grid points
    Param.lr = 101;
    Param.T = 1e-9;                   % Total simulation time [s]

    Param.dt = 1e-14;               % Time separation  [s]

    Param.V0 = 0;                    % Voltage at r=1 and t=1 [V]
    Param.VT = 0;                 % Ending voltage (r=end) [V]
    Param.EndV0 = -0.5;
    Param.EndVT = 0.8;


    Param.case = 5;
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian

end



