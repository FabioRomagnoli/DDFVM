function [alphalow, alphahigh] = alpha_exp(r, v, Ei, beta)
    Dv = diff(v);
    rhalf = (r(1:end-1) + r(2:end))/2;
    logPos = log(r(2:end) ./ r(1:end-1));

    eP = exp(-abs(Ei .* r(1:end-1).* logPos ./ (Dv-eps)));
    eN = exp(-abs(Ei .* r(2:end).* logPos ./ (Dv-eps)));
    eH = exp(-abs(Ei .* rhalf .* logPos ./ (Dv-eps)));
    
    alphalow  = beta * (eH - eN).*(Dv./(Ei*logPos));
    alphahigh = beta * (eP - eH).*(Dv./(Ei*logPos));
    
end
