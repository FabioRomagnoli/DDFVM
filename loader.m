clear all

loadSol = "nonConstFullResult.mat";
load(fullfile(".\sim\", loadSol));
fprintf("Loaded Solution %s\n",loadSol);



%%
file.Flag.concentrationPlot = "none";    % "all"/"last"/"none"
file.Flag.potentialPlot = "none";
file.Flag.currentPlot = "all";
file.Flag.generationPlot = "none";

plotter(file.Res,file.Dati,file.Flag);



%%


function concatenate(f1,f2)
    concat = f1.file;
    concat.Dati.T = concat.Dati.T + f2.file.Dati.T;
    concat.Dati.K = concat.Dati.K + f2.file.Dati.K;
    
    concat.Res.ASol = [concat.Res.ASol, f2.file.Res.ASol];
    concat.Res.Sol = [concat.Res.Sol, f2.file.Res.Sol];
    concat.Res.elapsedTime = concat.Res.elapsedTime + f2.file.Res.elapsedTime;
    concat.Res.kf = concat.Res.kf + f2.file.Res.kf;
    
    
    save(fullfile(".\sim\", "nonconstFullResult"), 'concat');
    fprintf("Saved Solution %a","nonconstFullResult");
end