function [F,jac] = assemblerPlasmaJGenAlphaExp(x, x0, BCs,  AD, dt)
    % Unpack
    M = AD.M;
    A = AD.A;
    A_bc = AD.A_bc;
    lrr = AD.lrr;
        
    % get boundary
    v_bc = BCs(1:2);
    n_bc = BCs(3:4);
    p_bc = BCs(5:6);

    % Extract vectors
    v = [v_bc(1); x(1:lrr); v_bc(2)];
    n = [n_bc(1); x(lrr+1:2*lrr); n_bc(2)];
    p = [p_bc(1); x(2*lrr+1:end); p_bc(2)];
    
    % Compute bimu bernulli only once
    DV =  diff(v) / AD.Vth;
    [Bp, Bn] = bimu_bernoulli (DV);

    % Matrix definitions
    An_full = ax_dd(AD.r, AD.mun, AD.Vth, Bn, Bp); % switched to account for z
    Ap_full = ax_dd(AD.r, AD.mup, AD.Vth, Bp, Bn);
   
    % Get reduced matrix
    An = An_full(2:end-1,2:end-1);
    Ap = Ap_full(2:end-1,2:end-1);
    
    % getting the first and last elements of the matrices 
    An_bc = An_full(2:end-1,[1 end]);
    Ap_bc = Ap_full(2:end-1,[1 end]);
    
    % Full matrix construction
    zeri = zeros(lrr);
    NL = [A M -M; 
        zeri (M+dt*An) zeri; 
        zeri zeri (M+dt*Ap)];
    
    % previous time step (no bounds)
    n0 = x0(lrr+1:2*lrr);
    p0 = x0(2*lrr+1:end);

    % system elements
    bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];

    rhs = [zeros(lrr,1);  M*n0  ; M*p0];

    % Generation Term
    Jn = comp_current(AD.r,AD.mun,AD.Vth,-1,n, Bn, Bp);
    [alphalow, alphahigh] = alpha_exp(AD.r, v, AD.Ei, AD.beta);
    R = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));
    gen = [zeros(lrr,1); dt*R(2:end-1) + dt*M*AD.S; dt*R(2:end-1) + dt*M*AD.S];

    
    % Build system
    F = NL*x - gen + bounds - rhs;

    % JACOBIAN 
    if nargout>1
        % JACOBIAN DOENS'T EXIST FOR THIS ASSEMBLER
    end
end
