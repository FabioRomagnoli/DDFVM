function [Param, Flag, Opt] = diodeCase3()

    % simulation settings
    Flag.model = "diode";
    Flag.verbose = 3;

    Flag.scheme = "coupled";
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)
    Flag.EndVT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)

    Flag.adaptive = false;        
    Flag.mesh = "linear";           % "linear"/"tanh"

    Flag.loadSol =  "no";
    Flag.saveSol = ".\diode\diodeCase4";
        
    
    % Hyperparameters Diode
    Param.K = 100;                   %Time grid points
    Param.lr = 101;
    Param.T = 0.0001;                   % Total simulation time [s]

    Param.dt = 1e-6;               % Time separation  [s]

    Param.V0 = -1.4;                    % Voltage at r=1 and t=1 [V]
    Param.VT = 1.4;                 % Ending voltage (r=end) [V]
    Param.EndV0 = 0;
    Param.EndVT = 0;


    Param.case = 3;
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian

end



