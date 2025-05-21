function [Param, Flag, Opt] = diodeCase3()

    % simulation settings
    Flag.model = "diode";
    Flag.verbose = 1;

    Flag.scheme = "coupled";
    Flag.method = "fsolve";         % "fsolve"/"newton"
    Flag.VT = "plateu";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)

    Flag.adaptive = false;        
    Flag.mesh = "linear";           % "linear"/"tanh"

    Flag.loadSol =  "no";
    Flag.saveSol = ".\diode\diodeReverse";
        
    

    % Hyperparameters Diode
    Param.K = 1000;                   %Time grid points
    Param.lr = 101;
    Param.T = 0.001;                   % Total simulation time [s]

    Param.dt = 1e-6;               % Time separation  [s]

    Param.V0 = 0;                    % Voltage at r=1 and t=1 [V]
    Param.VT = 0;                 % Ending voltage (r=end) [V]

   
    Param.EndV0 = 0;
    Param.EndVT = 1;


    Param.case = 2;
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian

end



