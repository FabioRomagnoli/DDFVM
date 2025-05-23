%% Plasma Characterstic Plot
clear all;
figure()
title('Applied Voltage vs measured Current')
hold on; 
    
dataPlourabu =  csvread("dati_esperimento_zheng\xx_PlourabuèData.csv");
dataPlourabu(all(dataPlourabu == 0, 2), :) = [];  % removes rows where both elements are zero
plot(dataPlourabu(:,1),dataPlourabu(:,2),"k-s", 'DisplayName', 'Experimental')

loadSol = "alphaExp50kV_51cells.mat";         % FULL 51 cells
load(fullfile(".\sim\", loadSol));
plot(file.Res.Vplot*1e-3,file.Res.JJ(1,:)*1e3, "r", 'linewidth',0.8,'DisplayName', '50 cells');

loadSol = "alphaExp50kV_81cells.mat";                 % FULL 81 cells
load(fullfile(".\sim\", loadSol));
plot(file.Res.Vplot*1e-3,file.Res.JJ(1,:)*1e3, 'Color', [0, 0.5, 0], 'linewidth',0.8, 'DisplayName', '80 cells');

loadSol = "alphaExp50kV_101cells.mat";           % FULL 101 cell
load(fullfile(".\sim\", loadSol));
plot(file.Res.Vplot*1e-3,file.Res.JJ(1,:)*1e3, "b", 'linewidth',0.8, 'DisplayName', '100 cells');
hold off;


grid on;
% set(gca, 'YScale', 'log') % Change y-axis to log scale
% set(gca, 'XScale', 'log') % Change y-axis to log scale

xlim([20 52]);
ylim([0 1.5]);

legend('Location', 'best'); 
ylabel("Current [mA/m]");
xlabel("Voltage [kV]");




%% Plasma concentration Plot
clear all;


loadSol = "alphaExp50kV_51cells.mat";         % FULL 51 cells
load(fullfile(".\sim\", loadSol));
vend = file.Res.Sol(file.Dati.vIdxs,end);

% 
% dataPositiveIon =  csvread("dati_esperimento_zheng\positive.csv");
% dataNodes =  csvread("dati_esperimento_zheng\nodes.csv");
% dataPotential =  csvread("dati_esperimento_zheng\potential.csv");
% 
% expVend = dataPotential(idx,:);
% 
% [~, idx] = min(abs(expVend(:,1) - Res.Sol(1,end)));
% 





