function plotPlasma(Sol, R, D, F)  
    % Concentration of each species
    concentrationPlot(Sol, R, D, F);

    % Potential V
    potentialPlot(Sol, R, D, F);

    % Electrical  current J
    currentPlot(Sol, R, D, F);

    % comparison of rhs generation, prev used for db
    generationPlot(Sol,R,D,F);

    % comparison of comp current methods, prev used for db
    compCurrentPlot(Sol,R,D,F)

    
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
    title('Current ')
    for k=ks:R.kf
        vk = Sol(D.vIdxs,k);
        nk = Sol(D.nIdxs,k);
        pk = Sol(D.pIdxs,k);

        clf; % Clear figure before plotting new data

        Jn = 2*pi*D.q*comp_current(D.r,D.mun,vk,D.Vth,-1,nk);
        Jp = 2*pi*D.q*comp_current(D.r,D.mup,vk,D.Vth, 1,pk);
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

function generationPlot(Sol, R, D, F)
    if F.generationPlot == "none", return; end

    % With matri
    Ar = ax_gen(D.r, Sol(D.vIdxs,end), D.mun, D.alpha, D.Vth, -1);
    
    % generation term
    dr = diff(D.r);
    Jn = comp_current(D.r,D.mun,Sol(D.vIdxs,end),D.Vth,-1,Sol(D.nIdxs,end));

    % Alpha term
    if strcmp(F.alpha,"const")
        alphalow = dr/2; alphahigh = dr/2;
    elseif strcmp(F.alpha,"exp")
        [alphalow, alphahigh] = alpha_exp(D.r, Sol(D.vIdxs,end), D.Ei, R.alpha);
    end
    
    RHS = abs(([0; alphalow.*Jn]+[alphahigh.*Jn; 0]));

    figure()
    title('Comparison of Integral of different generation terms')

    hold on; 
    
    plot(D.r,RHS,"b-o", 'DisplayName', 'Jn [1/m^2*s]');

    plot(D.r, R.genInt, "k-s", 'DisplayName', 'Const Gen');

    plot(D.r, Ar*Sol(D.nIdxs,end),"-*",'DisplayName', 'Matrix' );    
    hold off;

    grid on;
    set(gca, 'YScale', 'log') % Change y-axis to log scale
    set(gca, 'XScale', 'log') % Change y-axis to log scale
    xlim([D.r0 D.r0+D.ionLength]);
    legend('Location', 'best'); 
    ylabel("RHS [1/(ms)]");


end

function compCurrentPlot(Sol,R,D,F)
    if F.compCurrentPlot == "none", return; end

    rhalf = (D.r(1:end-1) + D.r(2:end))/2;
    dd = [0; diff(D.r)/2] + [diff(D.r)/2; 0];
    
    Arn = ax_gen(D.r, Sol(D.vIdxs,end), D.mun, 1, D.Vth, -1);
    AJn = ((2*pi*D.q)*Arn*Sol(D.nIdxs,end))./dd;

    Arp = ax_gen(D.r, Sol(D.vIdxs,end), D.mup, -1, D.Vth, -1);
    AJp = ((-2*pi*D.q)*Arp*Sol(D.pIdxs,end))./dd;

   
    figure()
    title('Comparison of Current computation')

    hold on; 
    plot(rhalf,R.Jn,"b-o", 'DisplayName', 'comp_current');
    plot(D.r,  AJn,"-*",'DisplayName', 'RHS' );    

    plot(rhalf,R.Jp,"b-o", 'DisplayName', 'comp_current');
    plot(D.r,  AJp,"-*",'DisplayName', 'RHS' );  
    hold off;

    grid on;
    set(gca, 'YScale', 'log') % Change y-axis to log scale
    set(gca, 'XScale', 'log') % Change y-axis to log scale
    legend('Location', 'best'); 

    % xlim([D.r0 D.r0+D.ionLength]);
    
end