function [Sol, Res] = postDiode(AD, D, ASol, F)
    % Redimensionalize solution 
    v = ASol(D.vIdxs,:)*AD.Vbar;
    n = ASol(D.nIdxs,:)*AD.nbar;
    p = ASol(D.pIdxs,:)*AD.nbar;
    Sol = [v; n; p];   % Sol will be 303 × 101
    
    Res.kf = length(Sol(1,:));
  
    Res.Jn = 2*pi*D.q*comp_current(D.r,D.mu,v(:,end),D.Vth,-1,n(:,end));
    Res.Jp = 2*pi*D.q*comp_current(D.r,D.mu,v(:,end),D.Vth, 1,p(:,end));
    Res.JJ = Res.Jn + Res.Jp;

    if std(Res.JJ) / mean(Res.JJ) < 1e-2 
        Ic = mean(Res.JJ);
        fprintf('Ic = %.5e', Ic);
    else
        fprintf('JJ is not constant. \n');
    end


    % Save file 
    if ~strcmp(F.saveSol,"no")
        saveFile.Sol = Sol;
        saveFile.ASol =  ASol;
        saveFile.Dati = Dati;
        saveFile.ADati = AD;
        saveFile.Flag = F;
        save(fullfile(".\sim\", F.saveSol), 'saveFile');
    end
end
