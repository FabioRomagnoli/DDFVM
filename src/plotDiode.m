loadSol = "diode\diodeSweep.mat";
% LOAD SOLUTION AND DEF PARAMETERS

load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s\n",loadSol);

% Get 4 linearly spaced indices from the full simulation range
numPlots = 4;
stepIndices = round(linspace(1, file.Res.kf, numPlots));

figure;
for i = 1:numPlots
    k = stepIndices(i);
    vk = file.Res.Sol(file.Dati.vIdxs, k);
    nk = file.Res.Sol(file.Dati.nIdxs, k);
    pk = file.Res.Sol(file.Dati.pIdxs, k);

    subplot(2,2,i);
    plot(file.Dati.r, pk, 'LineWidth', 1.2, 'DisplayName', 'p');
    hold on;
    plot(file.Dati.r, nk, 'LineWidth', 1.2, 'DisplayName', 'n');
    hold off;

    title(['V = ', num2str(file.Res.Vplot(1,k)), ' V']);
    xlabel('Position (m)');
    ylabel('Concentration (m^{-3})');
    legend('Interpreter', 'latex');
    % Uncomment the following to use log scale or limit axes
    % set(gca, 'YScale', 'log');
    % axis([Dati.r0, Dati.r1, 1e7, 1e15]);
    grid on;
end
