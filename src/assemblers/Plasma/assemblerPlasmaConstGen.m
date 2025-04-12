function [F,jac] = assemblerPlasmaConstGen(x, x0, BCs,  AD, dt)
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

    gen = [zeros(lrr,1); dt*M*AD.gen(2:end-1) ; dt*M*AD.gen(2:end-1) ];
    
    
    % Build system
    F = NL*x - gen + bounds - rhs;

    % JACOBIAN 
    if nargout>1
        zeri = zeros(lrr);
        
        dAn = zeros(AD.lr, AD.lr);    
        dAp = zeros(AD.lr, AD.lr);    
    
        % derivatives in dv
        for i=2:AD.lr-1
            log_pos = 1/log(AD.r(i)/AD.r(i+1));
            log_neg = 1/log(AD.r(i-1)/AD.r(i));
            
            up_neg = (v(i+1) - v(i))/AD.Vth;
            up_pos = (v(i) - v(i+1))/AD.Vth;
            dw_neg = (v(i) - v(i-1))/AD.Vth;
            dw_pos = (v(i-1) - v(i))/AD.Vth;
    
    
            dAn(i,i-1) = log_neg * (n(i)*DB(dw_neg) + n(i-1)*DB(dw_pos));
            dAn(i,i) = log_pos * (-n(i+1)*DB(up_neg) - n(i)*DB(up_pos)) +...
                            log_neg * (-n(i)*DB(dw_neg) - n(i-1)*DB(dw_pos));
            dAn(i,i+1) = log_pos * (n(i+1)*DB(up_neg) + n(i)*DB(up_pos));
    
    
            dAp(i,i-1) = log_neg * (-p(i)*DB(-dw_neg) - p(i-1)*DB(-dw_pos));
            dAp(i,i) = log_pos * (p(i+1)*DB(-up_neg) + p(i)*DB(-up_pos)) +...
                            log_neg * (p(i)*DB(-dw_neg) + p(i-1)*DB(-dw_pos));
            dAp(i,i+1) = log_pos * (-p(i+1) * DB(-up_neg) - p(i)*DB(-up_pos));
       
    
        end
    
        J11 = A;
        J12 = M;
        J13 = -M;
    
        J21 = dt*AD.mun*dAn(2:end-1, 2:end-1);
        J22 = M + dt*An;
        J23 = zeri;
    
        J31 = dt*AD.mup*dAp(2:end-1, 2:end-1); 
        J32 = zeri;
        J33 = M + dt*Ap;
    
        jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];    
    end
end
