function Ar = ax_gen(r, mu, alpha, Vth, z, Bp, Bn)
    lr = length(r);
    dr = diff(r);

    % Spdiags scarta primo elem sopradiag e ultimo sottodiag
    ld = [-dr/2.*Bp./log(r(1:end-1)./r(2:end)); 0];
    ud = [0; dr/2.*Bn./log(r(1:end-1)./r(2:end))];
    
    dd = ld+ud;

    % Dovrebbe avere un modulo, per questo ho aggiunto un -
    Ar = z*alpha.*(mu*Vth*spdiags([ld,dd, ud], -1:1, lr, lr));

end