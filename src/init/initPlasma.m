function  [Param,Dati, ADati] = initPlasma(Param, Flag)
    % Initialize missing parameters
    Param = paramPlasma(Param);
    
    % Initialize initial conditons
    Dati = datiPlasma(Param, Flag);
    
    % Adimensionalze all the parameters
    ADati = adimPlasma(Dati, Flag);
end


function Param = paramPlasma(Param)
    [def, Param] = plasmaDefaults(Param);

    % Get field names from the default struct
    defaultFields = fieldnames(def);
    % Loop over each field and add it to s_in if it's not already present
    for i = 1:length(defaultFields)
        if ~isfield(Param, defaultFields{i})
            Param.(defaultFields{i}) = def.(defaultFields{i});
        end
    end
end

function Dati = datiPlasma(Param, Flag)
    Dati = Param;

    % Mesh 
    delta   = 2.4577e-07; %2.4577e-05; % original value  
    alpha_msh  = -0.1;
    msh = CreateTanhMsh(Dati.lr, Dati.r0, Dati.r1 - Dati.r0, delta, alpha_msh);
    Dati.r = msh.x;
    Dati.lrr = Dati.lr - 2;

    % Vector
    Dati.S = Param.S*ones(Dati.lr-2,1);   
    
    % INITIAL CONDITIONS
    n0 = ones(Dati.lr,1) * Dati.N;
    p0 = ones(Dati.lr,1) * Dati.N;
    
    % solving for v0 inside domain 
    F = ax_laplacian (Dati.r, Dati.eps);    % Perché cambia il valore del residuo (quello restituito da fsolve) in base a se metto 1 o epsilon? l'rhs è 0!
    v0 = zeros(size(Dati.r));
    v0(1) = Dati.V0;               % can be inbetween 0 < v(0) < 5e4   [V]
    v0(end) = 0;
    % Siccome p-n+N fa 0 allora al rhs c'è solo la correzione con le condizioni al bordo
    v0(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * v0([1 end]) + Dati.q * (p0(2:Dati.lr-1) - n0(2:Dati.lr-1)));
    Dati.n0 = n0;  Dati.p0 = p0;  Dati.v0 = v0;
    
    Dati.x0 = [v0; n0; p0];
    % Internal index for full vectors
    Dati.vIdxs  = 1:Dati.lr;
    Dati.nIdxs  = Dati.lr + 1 : 2*Dati.lr;
    Dati.pIdxs  = 2*Dati.lr + 1 : 3*Dati.lr;

    % Internal indexes for reduced vectors 
    lrr = Dati.lr - 2;
    Dati.vIR  = 1:lrr;
    Dati.nIR  = lrr + 1 : 2*lrr;
    Dati.pIR  = 2*lrr + 1 : 3*lrr;

    % indexes for interior and boundary elements of full vector x 
    Dati.intIdxs = [2:Dati.lr-1, Dati.lr+2:2*Dati.lr-1, 2*Dati.lr+2:3*Dati.lr-1];
    Dati.bcsIdxs = [1, Dati.lr,  Dati.lr+1,  2*Dati.lr, 2*Dati.lr+1,  3*Dati.lr];

    if strcmp(Flag.VT, "linear")
        % linear increase from V0 to VT
        Dati.Vt = @(t) Dati.V0 + (Dati.VT - Dati.V0)/Dati.T*t;
    elseif strcmp(Flag.VT, "piecewise")
        % Attempt at piecewise V(t), steady for short period at V0 and VT
        Dati.Vt = @(t) (t <= Dati.T/100) * Dati.V0 + ...
                       (Dati.T/100 < t && t <= Dati.T-Dati.T/100) * (Dati.V0 + (Dati.VT - Dati.V0)/(Dati.T*0.98)*t) + ...
                       (Dati.T-Dati.T/100 < t )  * Dati.VT;
    elseif strcmp(Flag.VT, "plateu")
        Dati.Vt = @(t)  (t <= Dati.T-Dati.T/100) * (Dati.V0 + (Dati.VT - Dati.V0)/(Dati.T*0.99)*t) + ...
                        (Dati.T-Dati.T/100 < t )  * Dati.VT;
    end
    

    if strcmp(Flag.EndVT, "linear")
        % linear increase from V0 to VT
        Dati.EndVt = @(t) Dati.V0 + (Dati.EndVT - Dati.EndV0)/Dati.T*t;
    elseif strcmp(Flag.EndVT, "plateu")
        Dati.EndVt = @(t)  (t <= Dati.T-Dati.T/100) * (Dati.EndV0 + (Dati.EndVT - Dati.EndV0)/(Dati.T*0.99)*t) + ...
                        (Dati.T-Dati.T/100 < t )  * Dati.EndVT;
    end

    Dati.tsave = linspace(0, Dati.T, Dati.K+1);

    % Generation term
    Dati.ionLength = 1.4351e-04;  %1.4371e-04 is the true one but with lr=101 this is better
    idx = find(Dati.r >= Dati.r0 + Dati.ionLength,1);
    M_full = ax_mass(Dati.r, 1);
    genfull = zeros(Dati.lr,1);
    genfull(2:idx) = 1;
    G = (Dati.Iz/Dati.q)/(2*pi*sum(M_full*genfull));
    Dati.gen = G*genfull;
end

function AD = adimPlasma(D, Flag)
    AD = D;

    % Scaling factors
    AD.xbar = D.r1 - D.r0;                    % Lenght of device [m]
    AD.nbar = D.N;                            % not same meaning as diode [1/m^3]
    AD.Vbar = max(D.VT, D.Vth);                           % [V]
    AD.mubar= max(D.mun,D.mup);               % [m2]/([s][V])
    AD.tbar = 1/(AD.mubar*AD.Vbar/AD.xbar^2);          % [s]
    AD.Jbar = AD.mubar*AD.Vbar*AD.nbar/AD.xbar;           % [
    AD.qbar = D.q;


    % Adaptive time step parameters 
    AD.dtMin = D.dtMin/AD.tbar;
    AD.dtMax = D.dtMax/AD.tbar;

    % mesh
    AD.r1 = D.r1/AD.xbar; 
    AD.r0 = D.r0/AD.xbar;
    AD.r = D.r/AD.xbar;
    AD.ionLength = D.ionLength/AD.xbar;
    
    % time
    AD.T = D.T/AD.tbar;
    AD.dt = D.dt/AD.tbar; 
    AD.tsave = D.tsave/AD.tbar;

    % constants
    AD.N = D.N/AD.nbar;        % Useless
    AD.q = D.q/AD.qbar;        % maybe useless
    AD.Vth = D.Vth / AD.Vbar;
    AD.eps = D.eps*AD.Vbar/(D.q*AD.nbar*AD.xbar^2);     % adimensionale (squared normalized D. L.)
    AD.mun = D.mun/AD.mubar; 
    AD.mup = D.mup/AD.mubar;

    % plasma constants
    AD.S = D.S * AD.xbar^2/(AD.mubar*AD.Vbar*AD.nbar);
    AD.beta = D.beta * AD.xbar;
    AD.Ei = D.Ei * AD.xbar / AD.Vbar; 
    AD.alpha = D.alpha * AD.xbar;
    AD.gen = D.gen /(AD.Jbar/AD.xbar);

    % simulation parameters 
    AD.V0 = D.V0/AD.Vbar;
    AD.VT = D.VT/AD.Vbar;
    AD.Vz = D.Vz/AD.Vbar;
    AD.Iz = D.Iz*AD.xbar/AD.Vbar;

    if strcmp(Flag.VT, "linear")
        % linear increase from V0 to VT
        AD.Vt = @(t) AD.V0 +(AD.VT - AD.V0)/AD.T * t;
    elseif strcmp(Flag.VT, "piecewise")
        % Attempt at piecewise V(t), steady for short period at V0 and VT
        AD.Vt = @(t) (t <= AD.T/100) * AD.V0 + ...
                   (AD.T/100 < t && t <= AD.T-AD.T/100) * (AD.V0 + (AD.VT - AD.V0)/(AD.T*0.98)*t) + ...
                   (AD.T-AD.T/100 < t )  * AD.VT;
    elseif strcmp(Flag.VT, "plateu")
        AD.Vt = @(t)  (t <= AD.T-AD.T/100) * (AD.V0 + (AD.VT - AD.V0)/(AD.T*0.99)*t) + ...
                        (AD.T-AD.T/100 < t )  * AD.VT;
    end


    if strcmp(Flag.EndVT, "linear")
        % linear increase from V0 to VT
        AD.EndVt = @(t) AD.EndV0 +(AD.EndVT - AD.EndV0)/AD.T * t;
    elseif strcmp(Flag.EndVT, "plateu")
        AD.EndVt = @(t)  (t <= AD.T-AD.T/100) * (AD.EndV0 + (AD.EndVT - AD.EndV0)/(AD.T*0.99)*t) + ...
                        (AD.T-AD.T/100 < t )  * AD.EndVT;
    end

    % vectors
    AD.p0 = D.p0/AD.nbar;
    AD.n0 = D.n0/AD.nbar;
    AD.v0 = D.v0/AD.Vbar;
    AD.x0 = [AD.v0; AD.n0; AD.p0];
    
    % Concatenate all together into a single array
    AD.xBarVec = [repmat(AD.Vbar, 1, AD.lr), repmat(AD.nbar, 1, AD.lr), repmat(AD.nbar, 1, AD.lr)]';

    % Matrices
    AD.M_full = ax_mass(AD.r, 1);
    AD.M = AD.M_full(2:end-1,2:end-1);
    
    AD.A_full = ax_laplacian(AD.r,AD.eps);
    AD.A = AD.A_full(2:end-1,2:end-1);
    AD.A_bc = AD.A_full(2:end-1,[1 end]);
end


