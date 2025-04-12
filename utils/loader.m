clear all

loadSol = "alphaExpBeta7e5.full.mat";
load(fullfile(".\sim\", loadSol));
fprintf("Loaded Solution %s\n",loadSol);

%% 
file.Flag.saveSol = "no";
Res = postProcess(file.Dati, file.ADati, file.Res, file.Flag);

%%
file.Flag.concentrationPlot = "none";    % "all"/"last"/"none"
file.Flag.potentialPlot = "none";
file.Flag.currentPlot = "none";
file.Flag.generationPlot = "Y";

plotter(Res,file.Dati,file.Flag);


%%

% f1Name = "..\alphaExpBeta7e5.temp";
% f2Name = "alphaExpBeta7e5.3";
% saveName = "alphaExpBeta7e5.full.mat";

% concatenate(f1Name,f2Name,saveName);

function concatenate(f1Name,f2Name, saveName)
    f1 = load(fullfile(".\sim\pieces\", f1Name));
    f2 = load(fullfile(".\sim\pieces\", f2Name));

    file = f1.file;
    file.Dati.T = file.Dati.T + f2.file.Dati.T;
    % file.Dati.K = file.Dati.K + f2.file.Dati.K;
    
    file.Res.ASol = [file.Res.ASol, f2.file.Res.ASol];
    file.Res.Sol = [file.Res.Sol, f2.file.Res.Sol];
    file.Res.elapsedTime = file.Res.elapsedTime + f2.file.Res.elapsedTime;
    file.Res.kf = file.Res.kf + f2.file.Res.kf;
    
    
    save(fullfile(".\sim\", saveName), 'file');
    fprintf("Saved Solution %s\n",saveName);
end