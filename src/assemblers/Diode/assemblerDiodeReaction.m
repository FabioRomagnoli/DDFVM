function [F, jac]= assemblerDiodeReaction(x, x0, BCs,  AD, dt)
    % Unpack
    M = AD.M;
    lrr = AD.lrr;

    % previous time step (no bounds)
    n0 = x0(lrr+1:2*lrr);
    p0 = x0(2*lrr+1:end);

    R = (AD.ni^2-n0.*p0)./(AD.tau*(n0+p0));

    zeri = zeros(lrr);
    NL = [zeri  zeri  zeri; 
        zeri  M     zeri; 
        zeri  zeri  M];


    gen = [zeros(lrr,1); dt*M*R; dt*M*R];
    rhs = [zeros(lrr,1); M*n0; M*p0];
    
    F = NL*x - gen - rhs;

    if nargout>1
        J11 = zeri;
        J12 = zeri;
        J13 = zeri;
        J21 = zeri;
        J22 = M;
        J23 = zeri;
        J31 = zeri; 
        J32 = zeri;
        J33 = M;
        jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];

    end
end
