function Res = postProcess(D, AD, Res, Flag)
    Res = redim(D,AD,Res);

    % computes the final current if steady
    Res = finalCurrent(D, Res, Flag);

    if  strcmp(Flag.model, "plasma")
        Res = computeAlpha(D, Res, Flag);
    end

    saveFile(D, AD, Res, Flag)
end


function Res = redim(D,AD,Res)
    % Redimensionalize solution 
    v = Res.ASol(D.vIdxs,:)*AD.Vbar;
    n = Res.ASol(D.nIdxs,:)*AD.nbar;
    p = Res.ASol(D.pIdxs,:)*AD.nbar;
    Res.Sol = [v; n; p];   % Sol will be 303 × 101

    % number of elements in save file
    Res.kf = length(Res.Sol(1,:));
end


function Res = finalCurrent(D, Res, Flag)

    % Compute bimu bernulli only once
    DV =  diff(Res.Sol(D.vIdxs,end)) / D.Vth;
    [Bp, Bn] = bimu_bernoulli (DV);


    % the 2*pi comes from the radial contribution 
    Res.Jn = 2*pi*D.q*comp_current(D.r, D.mun, D.Vth, -1, Res.Sol(D.nIdxs,end), Bn, Bp);
    Res.Jp = 2*pi*D.q*comp_current(D.r, D.mup, D.Vth, 1, Res.Sol(D.pIdxs,end), Bp, Bn);
    Res.JJ = Res.Jn + Res.Jp;

    if std(Res.JJ) / mean(Res.JJ) < 1e-2 
        Ic = mean(Res.JJ);
        fprintf('\nIc = %.5e ', Ic);
        if strcmp(Flag.model,"plasma"), fprintf('Iz = %.5e, diff = %.5g',  D.Iz, abs(D.Iz - Ic)); end
    else
        fprintf('\nJJ is not constant.');
    end
end


function Res = computeAlpha(D, Res, Flag)
    if Flag.compAlpha
        Res.genInt = ax_mass(D.r, 1)*D.gen;
        intGen = sum(Res.genInt(2:end-1));

        % Compute bimu bernulli only once
        [Bp, Bn] = bimu_bernoulli (diff(Res.Sol(D.vIdxs,end)) / D.Vth);

        % generation term
        dr = diff(D.r);
        Jn = comp_current(D.r, D.mun, D.Vth, -1, Res.Sol(D.nIdxs,end), Bn, Bp);

        % Alpha term
        if strcmp(Flag.alpha,"const")
            alphalow = dr/2; alphahigh = dr/2;
        elseif strcmp(Flag.alpha,"exp")
            [alphalow, alphahigh] = alpha_exp(D.r, Res.Sol(D.vIdxs,end), D.Ei, 1);
        end
        
        R = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));
        intJ = sum(R(2:end-1));

        % Alpha is ratio of the two integrals over the ionization areas 
        if strcmp(Flag.alpha,"const")
            Res.alpha = intGen/intJ;
            Res.beta =  D.beta;
            fprintf('\nalpha = %.3f', Res.alpha);
        elseif strcmp(Flag.alpha,"exp")
            Res.beta =  intGen/intJ;
            Res.alpha = D.alpha;
            fprintf('\nbeta = %.3f', Res.beta);
        end
    else 
        Res.alpha = D.alpha;
        Res.beta = D.beta;
    end
end


function saveFile(D, AD, Res, Flag)
    if ~strcmp(Flag.saveSol,"no")
        file.Res = Res;
        file.Dati = D;
        file.ADati = AD;
        file.Flag = Flag;
        save(fullfile(".\sim\", Flag.saveSol), 'file');
        fprintf("\nSaved Solution to %s.m \n",Flag.saveSol);
    end
end