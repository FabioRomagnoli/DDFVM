function [alphalow, alphahigh] = alpha_exp(r, v, Ei, beta)
    Dv = diff(v);
    
    Dv(Dv == 0) = -eps;
    rhalf = (r(1:end-1) + r(2:end))/2;
    logPos = log(r(2:end) ./ r(1:end-1));

    eP = exp(Ei .* r(1:end-1).* logPos ./ Dv);
    eN = exp(Ei .* r(2:end).* logPos ./ Dv);
    eH = exp(Ei .* rhalf .* logPos ./ Dv);
    
    alphalow  = -beta * (eH - eN).*(Dv./(Ei*logPos));
    alphahigh = -beta * (eP - eH).*(Dv./(Ei*logPos));
end
