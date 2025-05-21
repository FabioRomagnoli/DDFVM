function plotter(Res, Dati, Flag)  
    % Concentration of each species
    concentrationPlot(Res, Dati, Flag);

    % Potential V
    potentialPlot(Res, Dati, Flag);

    % Electrical  current J
    currentPlot(Res, Dati, Flag);

    % comparison of rhs generation, prev used for debugging plasma
    generationPlot(Res,Dati,Flag);

    experimentalCurrentPotentialPlot(Res,Dati,Flag);
end


function concentrationPlot(Res, Dati, Flag)
    if isfield(Flag, 'concentrationPlot')
        switch Flag.concentrationPlot
            case "none"
                return
            case "last"
                ks = Res.kf;
                skip = 1;
            case "all"
                ks = 1;
                skip = 1;
            case "reduced"
                ks = 1;
                skip = 10;
        end    
    else
        return
    end

    figure;
    for k=ks:skip:Res.kf
        vk = Res.Sol(Dati.vIdxs,k);
        nk = Res.Sol(Dati.nIdxs,k);
        pk = Res.Sol(Dati.pIdxs,k);

        clf; % Clear figure before plotting new data
        title(['Electron and hole concentrations at V = ', num2str(vk(1))]);
        hold on;
        plot(Dati.r,pk, 'LineWidth', 1, "DisplayName", "p");
        plot(Dati.r,nk, 'LineWidth', 1, "DisplayName", "n");
        % plot(Dati.r,sqrt(nk.*pk), 'LineWidth',1, "DisplayName","$\sqrt{n+p}$");
        hold off;

        % grid on; 
        legend('Interpreter', 'latex');        
        % set(gca, 'YScale', 'log')
        % set(gca, 'XScale', 'log')
        % axis([Dati.r0,Dati.r1,1e7,1e15])
        xlabel('Position (m)');           % x-axis label
        ylabel('Concentration (m^{-3})'); % y-axis label


        drawnow;
    end
end

function potentialPlot(Res, Dati, Flag)
    if isfield(Flag, 'potentialPlot')
        switch Flag.potentialPlot
            case "none"
                return
            case "last"
                ks = Res.kf;
                skip = 1;
            case "all"
                ks = 1;
                skip = 1;
            case "reduced"
                ks = 1;
                skip = 10;
        end   
    else
        return
    end

    figure;
    title('Electric potential')
    for k = ks:skip:Res.kf

        clf; % Clear figure before plotting new data
        vk = Res.Sol(Dati.vIdxs, k);
        
        plot(Dati.r, vk, 'LineWidth', 1, "DisplayName", "v");
        
        grid on;
        legend();
        
        xlabel('Position (m)');         % x-axis label
        ylabel('Potential (V)');         % y-axis label
        
        drawnow;
    end
end
    
function currentPlot(Res, Dati, Flag)
    if isfield(Flag, 'currentPlot')
        switch Flag.currentPlot
            case "none"
                return
            case "last"
                ks = Res.kf;
                skip = 1;
            case "all"
                ks = 1;
                skip = 1;
            case "reduced"
                ks = 1;
                skip = 10;
            case "lastFew"
                ks = Res.kf - 100;
                skip = 1;
        end    
    else
        return
    end

    figure;

    rhalf = (Dati.r(1:end-1) + Dati.r(2:end))/2;
    for k=ks:skip:Res.kf

        vk = Res.Sol(Dati.vIdxs,k);
        nk = Res.Sol(Dati.nIdxs,k);
        pk = Res.Sol(Dati.pIdxs,k);

        % Compute bimu bernulli only once
        [Bp, Bn] = bimu_bernoulli (diff(vk) / Dati.Vth);

        clf; % Clear figure before plotting new data

        Jn = 2*pi*Dati.q*comp_current(Dati.r,Dati.mun,Dati.Vth,-1,nk, Bn, Bp);
        Jp = 2*pi*Dati.q*comp_current(Dati.r,Dati.mup,Dati.Vth, 1,pk, Bp, Bn);
        JJ = Jn + Jp;
        title(['Current at V = ', num2str(vk(1))]);

        hold on;
        plot(rhalf,Jn, 'LineWidth', 1, "DisplayName","Jn");
        plot(rhalf,Jp, 'LineWidth', 1, "DisplayName","Jp");
        plot(rhalf,JJ, 'LineWidth', 1, "DisplayName","JJ");
        hold off;

        grid on; 
        legend();
        % set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        drawnow;
    end
    
end

function generationPlot(Res, Dati, Flag)
    if ~isfield(Flag, 'generationPlot'), return; end
    if Flag.generationPlot == "none", return; end

    % Compute bimu_bernulli
    [Bp, Bn] = bimu_bernoulli (diff(Res.Sol(Dati.vIdxs,end)) / Dati.Vth);

    % With matrix
    Ar = ax_gen(Dati.r, Dati.mun, Res.alpha, Dati.Vth, -1, Bn, Bp);
    
    % generation term
    dr = diff(Dati.r);
    Jn = comp_current(Dati.r,Dati.mun,Dati.Vth,-1,Res.Sol(Dati.nIdxs,end),Bn, Bp);

    % Alpha term
    if strcmp(Flag.alpha,"const")
        alphalow = dr/2; alphahigh = dr/2;
    elseif strcmp(Flag.alpha,"exp")
        [alphalow, alphahigh] = alpha_exp(Dati.r, Res.Sol(Dati.vIdxs,end), Dati.Ei, Res.beta);
    end
    
    RHS = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));

    figure()
    title('Comparison of Integral of different generation terms')

    hold on; 
    
    plot(Dati.r, RHS,"b-o", 'DisplayName', 'Jn [1/m^2*s]');

    plot(Dati.r, Res.genInt, "k-s", 'DisplayName', 'Const Gen');

    plot(Dati.r, Ar*Res.Sol(Dati.nIdxs,end),"-*",'DisplayName', 'Matrix' );    
    hold off;

    grid on;
    set(gca, 'YScale', 'log') % Change y-axis to log scale
    set(gca, 'XScale', 'log') % Change y-axis to log scale
    % xlim([Dati.r0 Dati.r0+Dati.ionLength]);
    legend('Location', 'best'); 
    ylabel("RHS [1/(ms)]");

end



function experimentalCurrentPotentialPlot(Res, Dati, Flag)
    if ~isfield(Flag, 'experimentalCurrentPotentialPlot'), return; end
    if Flag.experimentalCurrentPotentialPlot == "none", return; end


    dataPlourabu =  csvread("dati_esperimento_zheng\xx_PlourabuèData.csv");
    dataPlourabu(all(dataPlourabu == 0, 2), :) = [];  % removes rows where both elements are zero

    figure()
    title('Applied Voltage vs measured Current plot')

    hold on; 

    plot(Res.Sol(1,:) - Res.Sol(Dati.lr,:),Res.JJv(:,1), "b", 'DisplayName', 'Simulation');

    plot(dataPlourabu(:,1)*1e3,dataPlourabu(:,2)*1e-3,"k-s", 'DisplayName', 'Experimental')

    hold off;

    grid on;
    % set(gca, 'YScale', 'log') % Change y-axis to log scale
    % set(gca, 'XScale', 'log') % Change y-axis to log scale

    % xlim([2e4 3.5e4]);
    % ylim([0 2e-4]);

    legend('Location', 'best'); 
    % ylabel("RHS [1/(ms)]");

end
