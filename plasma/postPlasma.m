function [Sol, Res] = postPlasma(AD, D, ASol, F)
    % Redimensionalize solution 
    v = ASol(D.vIdxs,:)*AD.Vbar;
    n = ASol(D.nIdxs,:)*AD.nbar;
    p = ASol(D.pIdxs,:)*AD.nbar;
    Sol = [v; n; p];   % Sol will be 303 × 101
    
    Res.kf = length(Sol(1,:));
  
    % the 2*pi comes from the radial contribution 
    Res.Jn = 2*pi*D.q*comp_current(D.r,D.mun,v(:,end),D.Vth,-1,n(:,end));
    Res.Jp = 2*pi*D.q*comp_current(D.r,D.mup,v(:,end),D.Vth, 1,p(:,end));
    Res.JJ = Res.Jn + Res.Jp;

    if std(Res.JJ) / mean(Res.JJ) < 1e-2 
        Ic = mean(Res.JJ);
        fprintf('Iz = %.5e, Ic = %.5e, diff = %.5g\n', D.Iz, Ic, abs(D.Iz - Ic));
    else
        fprintf('JJ is not constant. \n');
    end


    % Compute alpha 
    if strcmp(F.genterm, "const") && F.compAlpha
        Res.genInt = ax_mass(D.r, 1)*D.gen;
        intGen = sum(Res.genInt(2:end-1));
    
    
        % generation term
        dr = diff(D.r);
        Jn = comp_current(D.r,D.mun,v(:,end),D.Vth,-1,n(:,end));
    
        % Alpha term
        if strcmp(F.alpha,"const")
            alphalow = dr/2; alphahigh = dr/2;
        elseif strcmp(F.alpha,"exp")
            [alphalow, alphahigh] = alpha_exp(D.r, v(:,end), D.Ei, 1);
        end
        
        R = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));

        intJ = sum(R(2:end-1));

        % Alpha is ratio of the two integrals over the ionization areas 
        Res.alpha = intGen/intJ;
        fprintf('alpha = %.3f\n', Res.alpha);
    end


    % Save file 
    if ~strcmp(F.saveSol,"no")
        saveFile.Sol = Sol;
        saveFile.ASol =  ASol;
        saveFile.Dati = D;
        saveFile.ADati = AD;
        saveFile.Flag = F;
        save(fullfile(".\sim\", F.saveSol), 'saveFile');
        fprinf("Saved Solution %a",F.saveSol);

    end
end

