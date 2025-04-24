function [F,jac] = assemblerPlasmaReaction(x, x0, BCs,  AD, dt)
    % Unpack
    M = AD.M;
    lrr = AD.lrr;
        
    % previous time step (no bounds)
    n0 = x0(lrr+1:2*lrr);
    p0 = x0(2*lrr+1:end);

    % get boundary
    v_bc = BCs(1:2);
    n_bc = BCs(3:4);

    % Extract vectors
    v = [v_bc(1); x(1:lrr); v_bc(2)];
    n0full = [n_bc(1); n0; n_bc(2)];
    
    % Compute bimu bernulli only once
    DV =  diff(v) / AD.Vth;
    [Bp, Bn] = bimu_bernoulli (DV);

    % Full matrix construction
    zeri = zeros(lrr);
    NL = [zeri  zeri  zeri; 
        zeri  M     zeri; 
        zeri  zeri  M];
   

    rhs = [zeros(lrr,1);  M*n0  ; M*p0];

    % Generation Term
    Jn = comp_current(AD.r,AD.mun,AD.Vth,-1,n0full, Bn, Bp);
    [alphalow, alphahigh] = alpha_exp(AD.r, v, AD.Ei, AD.beta);
    R = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));
    gen = [zeros(lrr,1); dt*R(2:end-1) + dt*M*AD.S; dt*R(2:end-1) + dt*M*AD.S];

    
    % Build system
    F = NL*x - gen - rhs;

    % JACOBIAN 
    if nargout>1
        % JACOBIAN DOENS'T EXIST FOR THIS ASSEMBLER
    end
end
