function Res = postProcess(Dati, AD, Res, Flag, Param)
    Res = redim(Dati,AD,Res);

    if  strcmp(Flag.model, "plasma")
        Res = computeAlpha(Dati, Res, Flag);
    end

    if  Flag.computeCurrents
        Res = computeCurrents(Dati,Res, Flag, Param);
    end
    
    saveFile(Dati, AD, Res, Flag, Param)
end


function Res = redim(Dati,AD,Res)
    fprintf("\nRedimensionalized");

    % Redimensionalize solution 
    v = Res.ASol(Dati.vIdxs,:)*AD.Vbar;
    n = Res.ASol(Dati.nIdxs,:)*AD.nbar;
    p = Res.ASol(Dati.pIdxs,:)*AD.nbar;
    Res.Sol = [v; n; p];   % Sol will be 303 × 101

    % number of elements in save file
    Res.kf = length(Res.Sol(1,:));
end




function Res = computeCurrents(Dati,Res, Flag, Param)
    fprintf("\nComputed currents at every time step");

    for k=1:Res.kf
        vk = Res.Sol(Dati.vIdxs,k);
        nk = Res.Sol(Dati.nIdxs,k);
        pk = Res.Sol(Dati.pIdxs,k);
    
        % Compute bimu bernulli only once
        [Bp, Bn] = bimu_bernoulli (diff(vk) / Dati.Vth);

        % Compute current densities
        Jn = 2*pi * Dati.q * comp_current(Dati.r, Dati.mun, Dati.Vth, -1, nk, Bn, Bp);
        Jp = 2*pi * Dati.q * comp_current(Dati.r, Dati.mup, Dati.Vth,  1, pk, Bp, Bn);
    
        % Store individual and total current
        ResNew.Jn(:,k) = Jn;
        ResNew.Jp(:,k) = Jp;
        ResNew.JJ(k) = mean(Jn + Jp);
        ResNew.JJv(:,k) = Jn + Jp;
    end

    Res.Jn = ResNew.Jn;
    Res.Jp = ResNew.Jp;
    Res.JJ = ResNew.JJ;
    Res.JJv = ResNew.JJv;

    Res.DeviceLength = Dati.r1-Dati.r0;
    if strcmp(Flag.model, "diode")
        % The sign corrects for different orientation of the diode
        if Res.Sol(Dati.lr+1,1)  > Res.Sol(2*Dati.lr +1,1) % N pole on the inside 
            Res.signConv = -1; 
        else                                               % P pole on the inside
            Res.signConv = 1;
        end
    elseif strcmp(Flag.model, "plasma")
        Res.signConv = 1;
    end


    if std(Res.JJ(end)) / mean(Res.JJ(end)) < 1e-2 
        Ic = Res.DeviceLength*Res.signConv*mean(Res.JJ(end));
        fprintf('\nIc = %.5e ', Ic);
        if strcmp(Flag.model,"plasma"), fprintf('Iz = %.5e, diff = %.5g',  Dati.Iz, abs(Dati.Iz - Ic)); end
    else
        fprintf('\nJJ is not constant.');
    end


    Res.Vplot = (Param.Vbias -  Res.Sol(1,:) -  Res.Sol(Dati.lr,:));
    % device length to integrate, should be mass matrix
    Res.Iplot = Res.DeviceLength*Res.signConv*Res.JJ;

    if strcmp(Flag.model, "diode")
        % To get the bias voltage
        v = Res.Vplot(2:end);
        i = Res.Iplot(2:end);
        [~, idx] = min(abs(i));    % idx is the position in v/i where |I| is smallest
        % Return the corresponding voltage
        Res.Vbias = v(idx);
        fprintf('\nZero-current bias is approximately %.5g V\n', Res.Vbias);
    end
end




function Res = computeAlpha(Dati, Res, Flag)
    if Flag.compAlpha
        Res.genInt = ax_mass(Dati.r, 1)*Dati.gen;
        intGen = sum(Res.genInt(2:end-1));

        % Compute bimu bernulli only once
        [Bp, Bn] = bimu_bernoulli (diff(Res.Sol(Dati.vIdxs,end)) / Dati.Vth);

        % generation term
        dr = diff(Dati.r);
        Jn = comp_current(Dati.r, Dati.mun, Dati.Vth, -1, Res.Sol(Dati.nIdxs,end), Bn, Bp);

        % Alpha term
        if strcmp(Flag.alpha,"const")
            alphalow = dr/2; alphahigh = dr/2;
        elseif strcmp(Flag.alpha,"exp")
            [alphalow, alphahigh] = alpha_exp(Dati.r, Res.Sol(Dati.vIdxs,end), Dati.Ei, 1);
        end
        
        R = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));
        intJ = sum(R(2:end-1));

        % Alpha is ratio of the two integrals over the ionization areas 
        if strcmp(Flag.alpha,"const")
            Res.alpha = intGen/intJ;
            Res.beta =  Dati.beta;
            fprintf('\nalpha = %.3f', Res.alpha);
        elseif strcmp(Flag.alpha,"exp")
            Res.beta =  intGen/intJ;
            Res.alpha = Dati.alpha;
            fprintf('\nbeta = %.3f', Res.beta);
        end
    else 
        Res.alpha = Dati.alpha;
        Res.beta = Dati.beta;
    end
end


function saveFile(D, AD, Res, Flag, Param)
    if ~strcmp(Flag.saveSol,"no")
        file.Res = Res;
        file.Dati = D;
        file.ADati = AD;
        file.Flag = Flag;
        file.Param = Param;
        save(fullfile(".\sim\", Flag.saveSol), 'file');
        fprintf("\nSaved Solution to %s.m \n",Flag.saveSol);
    end
end