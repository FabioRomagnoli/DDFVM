function [F,jac] = assembler(x, x0, BCs,  AD, Flag, t, dt)
    % get boundary
    v_bc = BCs(1:2);
    n_bc = BCs(3:4);
    p_bc = BCs(5:6);

    % Extract vectors
    v = [v_bc(1); x(1:AD.lrr); v_bc(2)];
    n = [n_bc(1); x(AD.lrr+1:2*AD.lrr); n_bc(2)];
    p = [p_bc(1); x(2*AD.lrr+1:end); p_bc(2)];
    

    % Matrix definitions
    An_full = ax_dd(AD.r, v, AD.mun, AD.Vth, -1);
    Ap_full = ax_dd(AD.r, v, AD.mup, AD.Vth, 1); 
   
    % Get reduced matrix
    An = An_full(2:end-1,2:end-1);
    Ap = Ap_full(2:end-1,2:end-1);
    
    % getting the first and last elements of the matrices 
    A_bc = A_full(2:end-1,[1 end]);
    An_bc = An_full(2:end-1,[1 end]);
    Ap_bc = Ap_full(2:end-1,[1 end]);
    
    % Full matrix construction
    zeri = zeros(AD.lrr);
    NL = [A M -M; 
        zeri (M+dt*An) zeri; 
        zeri zeri (M+dt*Ap)];
    
    % previous time step (no bounds)
    n0 = x0(AD.lrr+1:2*AD.lrr);
    p0 = x0(2*AD.lrr+1:end);

    % system elements
    bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];

    if strcmp(Flag.model,"diode")
        rhs = [M*AD.N(2:end-1); M*n0; M*p0];
    elseif strcmp(Flag.model,"plasma")
        rhs = [zeros(AD.lrr,1);  M*n0  ; M*p0];
    end

    gen = generation(v, n, p, AD,Flag, dt, M);
    
    
    % Build system
    F = NL*x - gen + bounds - rhs;

    % JACOBIAN 
    if nargout>1
        jac = jacobian(v, n, p, AD, Flag, dt, A, M, An, Ap);
    end
end



function gen =  generation(v, n, p, AD, Flag, dt, M)
    if strcmp(Flag.model,"diode")
        R = @(n,p) (AD.ni^2 - n.*p)./(AD.tau*(n+p));
        gen = [zeros(AD.lrr,1); dt*M*R(n(2:end-1),p(2:end-1)); dt*M*R(n(2:end-1),p(2:end-1))];
   
    elseif strcmp(Flag.model,"plasma")
        if strcmp(Flag.genterm, 'non-const')
            % generation term
            dr = diff(AD.r);
            Jn = comp_current(AD.r,AD.mun,v,AD.Vth,-1,n);
        
            % Alpha term
            if strcmp(Flag.alpha,"const")
                alphalow = AD.alpha .* dr/2; alphahigh = AD.alpha .* dr/2;
            elseif strcmp(Flag.alpha,"exp")
                [alphalow, alphahigh] = alpha_exp(AD.r, v, AD.Ei, AD.beta);
            end
            
            R = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));
       
            gen = [zeros(AD.lrr,1); dt*R(2:end-1) + dt*M*AD.S; dt*R(2:end-1) + dt*M*AD.S];
        
        elseif strcmp(Flag.genterm, 'const')
            gen = [zeros(AD.lrr,1); dt*M*AD.gen(2:end-1) ; dt*M*AD.gen(2:end-1) ];
        end
        
    end
end


function jac = jacobian(v, n, p, AD, Flag, dt,  A, M, An, Ap)
    zeri = zeros(AD.lrr);
    
    dAn = zeros(AD.lr, AD.lr);    
    dAp = zeros(AD.lr, AD.lr);    
    dAr = zeros(AD.lr, AD.lr);    

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
   

        if strcmp(Flag.model,"plasma") && strcmp(Flag.genterm, 'non-const')
            dAr(i,i-1) = -((AD.r(i)-AD.r(i-1))/2) *log_neg * (n(i)*DB(dw_neg) + n(i-1)*DB(dw_pos));
            dAr(i,i) = ((AD.r(i+1)-AD.r(i))/2) * log_pos * (-n(i+1)*DB(up_neg) - n(i)*DB(up_pos)) -...
                            ((AD.r(i)-AD.r(i-1))/2) * log_neg * (-n(i)*DB(dw_neg) - n(i-1)*DB(dw_pos));
            dAr(i,i+1) = ((AD.r(i+1)-AD.r(i))/2) * log_pos * (n(i+1)*DB(up_neg) + n(i)*DB(up_pos));
        end
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


    if strcmp(Flag.model,"diode")
        dR_dn=@(n,p) -(p.^2 + AD.ni^2)./(AD.tau*((n+p).^2));
        dR_dp=@(n,p) -(n.^2 + AD.ni^2)./(AD.tau*((n+p).^2));
    
        dRn = diag(dt*M*dR_dn(n(2:end-1),p(2:end-1)));  
        dRp = diag(dt*M*dR_dp(n(2:end-1),p(2:end-1))); 

        J22 = J22 - dRn;
        J23 = J23 - dRp;

        J32 = J32 - dRn;
        J33 = J33 - dRp;

    elseif strcmp(Flag.model,"plasma") && strcmp(Flag.genterm, 'non-const')
        Ar_full = ax_gen(AD.r, v, AD.mun, AD.alpha, AD.Vth, -1); % generation term
        dRn = dt*Ar_full(2:end-1,2:end-1);

        J21 = J21 - AD.alpha*AD.mun*dt*dAr(2:end-1,2:end-1);
        J22 = J22 - dRn;
        
        J31 = J31 - AD.alpha*AD.mun*dt*dAr(2:end-1,2:end-1); 
        J32 = J32 - dRn;
    end


    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];
end 