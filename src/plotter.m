function plotter(R, D, F)  
    % Concentration of each species
    concentrationPlot(R, D, F);

    % Potential V
    potentialPlot(R, D, F);

    % Electrical  current J
    currentPlot(R, D, F);

    % comparison of rhs generation, prev used for debugging plasma
    generationPlot(R,D,F);
end


function concentrationPlot(R, D, F)
    if isfield(F, 'concentrationPlot')
        switch F.concentrationPlot
            case "none"
                return
            case "last"
                ks = R.kf;
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
    for k=ks:skip:R.kf
        vk = R.Sol(D.vIdxs,k);
        nk = R.Sol(D.nIdxs,k);
        pk = R.Sol(D.pIdxs,k);

        clf; % Clear figure before plotting new data
        title(['Electron and hole concentrations at V = ', num2str(vk(1))]);
        hold on;
        plot(D.r,pk, 'LineWidth', 1, "DisplayName", "p");
        plot(D.r,nk, 'LineWidth', 1, "DisplayName", "n");
        % plot(D.r,sqrt(nk.*pk), 'LineWidth',1, "DisplayName","$\sqrt{n+p}$");
        hold off;

        % grid on; 
        legend('Interpreter', 'latex');        
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')

        xlabel('Position (m)');           % x-axis label
        ylabel('Concentration (m^{-3})'); % y-axis label


        drawnow;
    end
end

function potentialPlot(R, D, F)
    if isfield(F, 'potentialPlot')
        switch F.potentialPlot
            case "none"
                return
            case "last"
                ks = R.kf;
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
    for k = ks:skip:R.kf

        clf; % Clear figure before plotting new data
        vk = R.Sol(D.vIdxs, k);
        
        plot(D.r, vk, 'LineWidth', 1, "DisplayName", "v");
        
        grid on;
        legend();
        
        xlabel('Position (m)');         % x-axis label
        ylabel('Potential (V)');         % y-axis label
        
        drawnow;
    end
end
    
function currentPlot(R, D, F)
    if isfield(F, 'currentPlot')
        switch F.currentPlot
            case "none"
                return
            case "last"
                ks = R.kf;
                skip = 1;
            case "all"
                ks = 1;
                skip = 1;
            case "reduced"
                ks = 1;
                skip = 10;
            case "lastFew"
                ks = R.kf - 100;
                skip = 1;
        end    
    else
        return
    end

    figure;

    rhalf = (D.r(1:end-1) + D.r(2:end))/2;
    for k=ks:skip:R.kf

        vk = R.Sol(D.vIdxs,k);
        nk = R.Sol(D.nIdxs,k);
        pk = R.Sol(D.pIdxs,k);

        % Compute bimu bernulli only once
        [Bp, Bn] = bimu_bernoulli (diff(vk) / D.Vth);

        clf; % Clear figure before plotting new data
        title(['Current at V = ', num2str(vk(1))]);

        Jn = 2*pi*D.q*comp_current(D.r,D.mun,D.Vth,-1,nk, Bn, Bp);
        Jp = 2*pi*D.q*comp_current(D.r,D.mup,D.Vth, 1,pk, Bp, Bn);
        JJ = Jn + Jp;

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

function generationPlot(R, D, F)
    if ~isfield(F, 'generationPlot'), return; end
    if F.generationPlot == "none", return; end

    % Compute bimu_bernulli
    [Bp, Bn] = bimu_bernoulli (diff(R.Sol(D.vIdxs,end)) / D.Vth);

    % With matrix
    Ar = ax_gen(D.r, D.mun, R.alpha, D.Vth, -1, Bn, Bp);
    
    % generation term
    dr = diff(D.r);
    Jn = comp_current(D.r,D.mun,D.Vth,-1,R.Sol(D.nIdxs,end),Bn, Bp);

    % Alpha term
    if strcmp(F.alpha,"const")
        alphalow = dr/2; alphahigh = dr/2;
    elseif strcmp(F.alpha,"exp")
        [alphalow, alphahigh] = alpha_exp(D.r, R.Sol(D.vIdxs,end), D.Ei, R.beta);
    end
    
    RHS = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));

    figure()
    title('Comparison of Integral of different generation terms')

    hold on; 
    
    plot(D.r, RHS,"b-o", 'DisplayName', 'Jn [1/m^2*s]');

    plot(D.r, R.genInt, "k-s", 'DisplayName', 'Const Gen');

    plot(D.r, Ar*R.Sol(D.nIdxs,end),"-*",'DisplayName', 'Matrix' );    
    hold off;

    grid on;
    set(gca, 'YScale', 'log') % Change y-axis to log scale
    set(gca, 'XScale', 'log') % Change y-axis to log scale
    xlim([D.r0 D.r0+D.ionLength]);
    legend('Location', 'best'); 
    ylabel("RHS [1/(ms)]");

end
