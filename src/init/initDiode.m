function  [Param,Dati, ADati] = initDiode(Param, Flag)
    % Initialize missing parameters
    Param = paramDiode(Param);
    
    % Initialize initial conditons
    Dati = datiDiode(Param, Flag);
    
    % Adimensionalze all the parameters
    ADati = adimDiode(Dati, Flag);
end

function Param = paramDiode(Param)
    [def, Param] = diodeDefaults(Param);

    % Get default values into Param
    defaultFields = fieldnames(def);
    % Loop over each field and add it to s_in if it's not already present
    for i = 1:length(defaultFields)
        if ~isfield(Param, defaultFields{i})
            Param.(defaultFields{i}) = def.(defaultFields{i});
        end
    end
end

function Dati = datiDiode(Param, Flag)
    Dati = Param;

    % Mesh
    if strcmp(Flag.mesh, "linear")
        Dati.r = linspace(Dati.r0,Dati.r1,Dati.lr)';   % grid points       [m]
    elseif strcmp(Flag.mesh, "tanh")
        % delta   = 2.4577e-05; % original value  
        delta   = 2.4577e-09;
        alpha_msh  = -0.1;
        msh = CreateTanhMsh(Dati.lr, Dati.r0, Dati.r1 - Dati.r0, delta, alpha_msh);
        Dati.r = msh.x;
    end
    Dati.lrr = Dati.lr - 2;

    % Vector
    Dati.N = Param.N*ones(size(Dati.r));        % N density of the gas      [m-3]
    
    % INITIAL CONDITIONS
    if Param.case == 1
        n0 = (Dati.N + sqrt(Dati.N.^2 + Dati.ni^2*4))/2;   
        p0 = Dati.ni^2./n0; 

    elseif Param.case == 2
        Nd  = 1e24;
        Na = 1e22;

        Dati.N = Nd * (Dati.r <= 1 + 500e-9) -...
            Na * (Dati.r > 1 + 500e-9);

        n0 = (Dati.N + sqrt(Dati.N.^2 + Dati.ni^2*4))/2;% Valore esatto dell'n (si trova dal fatto che p=ni^2/n
        p0 = (-Dati.N + sqrt(Dati.N.^2 + 4*Dati.ni^2))/2;
        p0(Dati.N > 0) = Dati.ni^2 ./ n0(Dati.N > 0); % Capire perché è necessario, forse è solo un problema di approx numerica, infatti N^2 è molto maggiore di ni^2
        n0(Dati.N < 0) = Dati.ni^2 ./ p0(Dati.N < 0);

    elseif Param.case == 3
        Nd  = 1e24;
        Na1 = 1e22;
        Na2 = 1e24;
        
        Dati.N = Nd * (Dati.r <= 1+50e-9) ...
            - Na1 * (Dati.r > 1+50e-9 & Dati.r <=1+ 500e-9) ...
            - Na2 * (Dati.r > 1+500e-9);
        
        n0 = (Dati.N + sqrt(Dati.N.^2 + Dati.ni^2*4))/2; % Valore esatto dell'n (si trova dal fatto che p=ni^2/n
        p0 = (-Dati.N + sqrt(Dati.N.^2 + 4*Dati.ni^2))/2;
        p0(Dati.N > 0) = Dati.ni^2 ./ n0(Dati.N > 0); % Capire perché è necessario, forse è solo un problema di approx numerica, infatti N^2 è molto maggiore di ni^2
        n0(Dati.N < 0) = Dati.ni^2 ./ p0(Dati.N < 0);
    elseif Param.case == 4
        Nd  = 1e24;
        Na = 1e24;

        Dati.N = Nd * (Dati.r <= 1 + 500e-9) -...
            Na * (Dati.r > 1 + 500e-9);

        n0 = (Dati.N + sqrt(Dati.N.^2 + Dati.ni^2*4))/2;% Valore esatto dell'n (si trova dal fatto che p=ni^2/n
        p0 = (-Dati.N + sqrt(Dati.N.^2 + 4*Dati.ni^2))/2;
        p0(Dati.N > 0) = Dati.ni^2 ./ n0(Dati.N > 0); % Capire perché è necessario, forse è solo un problema di approx numerica, infatti N^2 è molto maggiore di ni^2
        n0(Dati.N < 0) = Dati.ni^2 ./ p0(Dati.N < 0);
    
    elseif Param.case == 5
        Na = 1e24; % [1/m^3] - P-side (inner part)
        Nd = 1e22; % [1/m^3] - N-side (outer part)
    
        % Define junction at r = 1 + 500e-9 m
        % P-side for r <= 1 + 500e-9
        % N-side for r >  1 + 500e-9
        Dati.N = Nd * (Dati.r > 1 + 500e-9) - Na * (Dati.r <= 1 + 500e-9);
    
        % Equilibrium carrier concentrations
        n0 = (Dati.N + sqrt(Dati.N.^2 + 4 * Dati.ni^2)) / 2;
        p0 = (-Dati.N + sqrt(Dati.N.^2 + 4 * Dati.ni^2)) / 2;
    
        % Ensure numerical stability in high doping regions
        p0(Dati.N > 0) = Dati.ni^2 ./ n0(Dati.N > 0); % N-side
        n0(Dati.N < 0) = Dati.ni^2 ./ p0(Dati.N < 0); % P-side
    elseif Param.case == 6
        %....
    end




    % solving for v0 inside domain 
    F = ax_laplacian (Dati.r, 1);    % Perché cambia il valore del residuo (quello restituito da fsolve) in base a se metto 1 o epsilon? l'rhs è 0!
    v0 = zeros(size(Dati.r));
    v0(1) = Param.V0;               % can be inbetween 0 < v(0) < 5e4   [V]
    v0(end) = Param.EndV0;
    % Siccome p-n+N fa 0 allora al rhs c'è solo la correzione con le condizioni al bordo
    v0(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * v0([1 end]));
    Dati.n0 = n0;  Dati.p0 = p0;  Dati.v0 = v0;
    
    Dati.x0  = [v0; n0; p0];
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

end


function AD = adimDiode(D, Flag)
    AD = D;

    % Scaling factors
    AD.xbar = D.r1 - D.r0;                    % Lenght of device [m]
    AD.nbar = norm(D.N,'inf');                            % not same meaning as diode [1/m^3]
    AD.Vbar = D.Vth;                           % [V]
    AD.tbar = D.tau;          % [s]
    AD.mubar= AD.tbar*AD.Vbar/AD.xbar^2;               % [m2]/([s][V])
    AD.qbar = D.q;
    

    % Adaptive time step parameters 
    AD.dtMin = D.dtMin/AD.tbar;
    AD.dtMax = D.dtMax/AD.tbar;

    % mesh
    AD.r1 = D.r1/AD.xbar; 
    AD.r0 = D.r0/AD.xbar;
    AD.r = D.r/AD.xbar;
    
    % time
    AD.T = D.T/AD.tbar;
    AD.dt = D.dt/AD.tbar; 
    AD.tsave = D.tsave/AD.tbar;
    AD.tau = D.tau/AD.tbar;
    
    % constants
    AD.N = D.N/AD.nbar;        % Useless
    AD.q = D.q/AD.qbar;        % maybe useless
    AD.Vth = D.Vth / AD.Vbar;
    AD.eps = D.eps*AD.Vbar/(D.q*AD.nbar*AD.xbar^2);     % adimensionale (squared normalized D. L.)
    AD.ni = D.ni/AD.nbar;
    AD.mu = D.mu*AD.mubar; 
    AD.mun =  D.mun*AD.mubar;
    AD.mup = AD.mup*AD.mubar;

    % simulation parameters 
    AD.V0 = D.V0/AD.Vbar;
    AD.VT = D.VT/AD.Vbar;
    AD.EndV0 = D.EndV0/AD.Vbar;
    AD.EndVT = D.EndVT/AD.Vbar;


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
    AD.v0 = D.v0/AD.Vbar;
    AD.n0 = D.n0/AD.nbar;
    AD.p0 = D.p0/AD.nbar;
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