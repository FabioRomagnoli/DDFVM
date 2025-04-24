function [Param, Flag, Opt] = loader()

    Flag.loadSol = "plasmaConstGen.mat";
    Flag.saveSol = "no";

    Flag.init = false;
    Flag.solve = false;
    Flag.postProcess = false;
    
    Flag.concentrationPlot = "last";    % "all"/"last"/"none"
    Flag.potentialPlot = "none";
    Flag.currentPlot = "reduced";
    Flag.generationPlot = "none";

    Param = struct();
    Opt = struct();
end