function [def, Param] = plasmaDefaults(Param)    
    def = struct();

    % Time discretization
    if ~isfield(Param,"K"), Param.("K") = 100; end
    if ~isfield(Param,"dt"), Param.("dt") = 1e-4; end
    def.T = Param.K * Param.dt;        % Simulation time               
    
    % Load exp data from Zheng
    if ~isfield(Param,"ZhengIdx"), Param.("ZhengIdx") = 14; end
    load('.\utils\ExpData.mat');
    def.Iz = xxZheng(Param.ZhengIdx,"I").(1) * 10^(-4);     % microV/cm converted to V/m
    def.Vz = xxZheng(Param.ZhengIdx,"V").(1) * 10^(3);      % kV converted to V
    def.V0 = def.Vz;                   % Voltage at r=1 and t=1  [V]
    def.VT = def.Vz;   
    
    % PHYSICAL PARAMETERS 
    def.eps = 8.8e-12;             % Permittivity              [C]/([V][m])
    def.mup = 2e-4;                    % Diffusivity coefficients  [m2]/([s][V])
    def.mun = 1e2 * def.mup; 
    def.q = 1.6e-19;                   % Charge                    [C] 
    def.Vth = 26e-3;                   % Thermal voltage           [V]
    def.N = 1e7;                      % density constant [m-3]
    
    % Plasma related constants
    def.alpha = 11586.738;  
    def.beta = 7.2e5;         % impact ionization coefficient                    [m]/([s][V])
    def.Ei = 2.09e7;  
    def.S = 1e9;                      % random constant  [?]
    
    % Mesh
    def.lr = 101;
    def.r0 = 700e-6;
    def.r1 = def.r0 + 10.35e-2;
    
    % Newton Parameters
    def.Nmaxit = 150; 
    def.Ntoll = 1e-5;
    
    % Adaptive time step parameters
    def.tolError = 1e-3;
    def.stepBuffer = 0.9;
    def.scalingPower = 0.5;
    def.dtMin = 1e-15;
    def.dtMax = 1e-4;
end