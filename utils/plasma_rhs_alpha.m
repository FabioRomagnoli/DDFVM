function R = plasma_rhs_alpha(r, v,n , mu, Vth, z, Ei, beta)

    dr = diff(r);
    Dv = diff(v);
    rhalf = r(1:end-1) + dr;

    [Bp, Bn] = bimu_bernoulli (z .* Dv / Vth);
 
    nBp = n(1:end-1).*Bp;
    nBn = n(2:end).*Bn;
    
    logUp = log(r(1:end-1) ./ r(2:end));
    logDw = log(r(2:end) ./ r(1:end-1));
    logDw = logUp;
    
    eP = exp(-Ei .* r(1:end-1).* logDw ./ Dv);
    eN = exp(-Ei .* r(2:end).* logDw ./ Dv);
    eH = exp(-Ei .* rhalf .* logDw ./ Dv);

    xJ = mu*Vth ./ logUp .*(nBn-nBp);
    
    rhs_l = xJ .* (eH - eN).*(Dv./(Ei*logDw));
    rhs_u = xJ .* (eP - eH).*(Dv./(Ei*logDw));
    
    R = -beta*([0; rhs_l]+[rhs_u; 0]);

end