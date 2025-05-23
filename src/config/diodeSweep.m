function [Param, Flag, Opt] = diodeCase3()

    % simulation settings
    Flag.model = "diode";
    Flag.verbose = 4;

    Flag.scheme = "coupled";
    Flag.method = "fsolve";         % "fsolve"/"newton"
    
    Flag.VT = "piecewise";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)
    Flag.EndVT = "piecewise";                 % "linear"/"piecewise"/"plateu"  controls how VT changes in time    (init.m)

    Flag.adaptive = false;        
    Flag.mesh = "linear";           % "linear"/"tanh"


    Flag.loadSol =  "no";
    Flag.saveSol = ".\diode\diodeSweep";
        
    
    % Hyperparameters Diode
    Param.K = 50;                   %Time grid points
    Param.lr = 101;
    Param.T = 1e-3;                   % Total simulation time [s]
    Param.dt = 1e-5;               % Time separation  [s]

    % Param.dtMin = 1e-15;
    Param.dtMax = 1e-5;

    
    % TO FIND THE BIAS do a forward bias run. (
    % Param.case = 5;               % higher p on the inside
    % Param.Vbias = -0.86869;
    % Param.V0 = Param.Vbias; 
    % Param.VT = Param.Vbias;
    % Param.EndV0 = 0;              % N pole 
    % Param.EndVT = -1;              % N pole en
    % % Decrease EdnVT forc  forward bias, 
    % % inrease  it for reverse

   
    % Param.case = 2;                 % higher n on the inside
    % Param.Vbias = -0.86352;
    % Param.V0 = 0;              % N pole 
    % Param.VT = -1;              % N pole en
    % Param.EndV0 = Param.Vbias; 
    % Param.EndVT = Param.Vbias;
    % Dicrease VT for forwards
    % incease for reverse


    % Param.case = 3;                 % higher n on the inside
    % Param.Vbias = -1;
    % Param.V0 = 0;              % N pole 
    % Param.VT = -1;              % N pole en
    % Param.EndV0 = Param.Vbias; 
    % Param.EndVT = Param.Vbias;
    % Dicrease VT for forwards
    % incease for reverse


    Param.case = 4;                 % higher n on the inside
    Param.Vbias = -0.73673;
    % Param.Vbias = 0;

    Param.V0 = 0;              % N pole 
    Param.VT = -1;              % N pole en
    Param.EndV0 = Param.Vbias; 
    Param.EndVT = Param.Vbias;
    % Dicrease VT for forwards
    % incease for reverse
    
    % Fsolve flags
    Opt.Display = "off";                    % "off"/"iter"/"final"/"final-detailed"     
    Opt.SpecifyObjectiveGradient = true;    % Jacobian or no Jacobian
    Opt.FunctionTolerance = 1e-4;
end



