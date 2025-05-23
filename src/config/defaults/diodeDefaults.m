function [def, Param] = diodeDefaults(Param)
    def = struct();
    % Time Discretization
    if ~isfield(Param,"K"), Param.("K") = 100; end
    if ~isfield(Param,"dt"), Param.("dt") = 1e-10; end
    def.T = Param.K * Param.dt;        % Simulation time        
    
    % PHYSICAL PARAMETERS 
    def.ni = 6.14e+15;                    % Intrinsic concentration
    def.eps = 1.035e-10;               % Permittivity              [C]/([V][m])
    def.q = 1.6e-19;                   % Charge                    [C] 
    def.mu = 0.1;                    % Diffusivity coefficients  [m2]/([s][V])
    def.mun = 0.135;
    def.mup = 0.048;
    def.Vth = 26e-3;   
    def.tau = 1e-11;
    
    def.V0 = 1.4;                   % Starting voltage (r=1) [V]
    def.VT = 1.4;                   % Ending voltage (r=end) [V]
    def.EndV0 = 0;
    def.EndVT = 0;


    % Diode related constants
    def.N = 10^22;                    % density constant [m-3]
    def.case = 3;
    
    % Mesh
    def.lr = 101;
    def.r0 = 1;
    def.r1 = def.r0 + 1e-6;
    
    % Newton Parameters
    def.Nmaxit = 150; 
    def.Ntoll = 1e-5;
    
    % Adaptive time step parameters
    def.tolError = 1e-3; 
    def.stepBuffer = 0.8;
    def.scalingPower = 0.3;
    def.dtMin = 1e-17;
    def.dtMax = 1e-2;
end