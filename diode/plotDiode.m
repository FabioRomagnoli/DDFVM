function plotDiode(Sol, R, D, F)  
    % Concentration of each species
    concentrationPlot(Sol, R, D, F);

    % Potential  V
    potentialPlot(Sol, R, D, F);

    % Electrical  current J
    currentPlot(Sol, R, D, F);
end


function concentrationPlot(Sol, R, D, F)
    switch F.concentrationPlot
        case "none"
            return
        case "last"
            ks = R.kf;
        case "all"
            ks = 1;
    end
    figure;
    title('Electron and hole concentrations')
    for k=ks:R.kf
        nk = Sol(D.nIdxs,k);
        pk = Sol(D.pIdxs,k);

        clf; % Clear figure before plotting new data
        hold on;
        plot(D.r,pk, 'LineWidth', 1, "DisplayName", "p");
        plot(D.r,nk, 'LineWidth', 1, "DisplayName", "n");
        plot(D.r,sqrt(nk.*pk),'LineWidth', 1, "DisplayName", 'sqrt(n*p)' )
        hold off;

        grid on; 
        legend();
        set(gca, 'YScale', 'log')
        % set(gca, 'XScale', 'log')

        xlabel('Position (m)');           % x-axis label
        ylabel('Concentration (m^{-3})'); % y-axis label


        drawnow;
    end
end

function potentialPlot(Sol, R, D, F)
    switch F.potentialPlot
        case "none"
            return
        case "last"
            ks = R.kf;
        case "all"
            ks = 1;
    end        
    figure;
    title('Electric potential')
    for k = ks:R.kf
        clf; % Clear figure before plotting new data
        vk = Sol(D.vIdxs, k);
        
        plot(D.r, vk, 'LineWidth', 1, "DisplayName", "v");
        
        grid on;
        legend();
        
        xlabel('Position (m)');         % x-axis label
        ylabel('Potential (V)');         % y-axis label
        
        drawnow;
    end
end
    
function currentPlot(Sol, R, D, F)
    switch F.currentPlot
        case "none"
            return
        case "last"
            ks = R.kf;
        case "all"
            ks = 1;
    end    

    figure;

    rhalf = (D.r(1:end-1) + D.r(2:end))/2;
    title('Current')
    for k=ks:R.kf
        vk = Sol(D.vIdxs,k);
        nk = Sol(D.nIdxs,k);
        pk = Sol(D.pIdxs,k);

        clf; % Clear figure before plotting new data

        Jn = 2*pi*D.q*comp_current(D.r,D.mu,vk,D.Vth,-1,nk);
        Jp = 2*pi*D.q*comp_current(D.r,D.mu,vk,D.Vth, 1,pk);
        JJ = Jn + Jp;

        hold on;
        plot(rhalf,Jn, 'LineWidth', 1, "DisplayName","Jn");
        plot(rhalf,Jp, 'LineWidth', 1, "DisplayName","Jp");
        plot(rhalf,JJ, 'LineWidth', 1, "DisplayName","JJ");
        hold off;

        grid on; 
        legend();
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        drawnow;
    end
    
end
