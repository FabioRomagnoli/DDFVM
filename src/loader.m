% clear all

% PLASMA CASES
% FULL SIMULATIONS 
loadSol = "alphaExp50kV_51cells.mat";         % FULL 51 cells
% loadSol = "alphaExp50kV_81cells.mat";                 % FULL 81 cells
% loadSol = "alphaExp50kV_101cells.mat";           % FULL 101 cell

% DIODE CASES
% loadSol = "diode\diodeSweepCase4Forward";
% loadSol = "diode\diodeSweepCase4ForwardSplit";
% loadSol = "diode\diodeSweepCase4ForwardCoupled";
% loadSol = "diode\diodeSweepCase4ForwardSplitSub";


% loadSol = "diode\diodeSweepCase4Reverse";


% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s\n",loadSol);

% OUTPUT CONFIGURATON
printConfiguration(struct(),file.Param,struct(),file.Flag,struct(),struct());
fprintf("'\nElapsed Time: %g", file.Res.elapsedTime);

% % POST PROCESSING IF NECESSARY
file.Flag.saveSol = loadSol;     % !!Risk of overriding   /   loadSol
file.Flag.computeCurrents = true;
Res = postProcess(file.Dati, file.ADati, file.Res, file.Flag, file.Param);


% PLOTTING
% file.Flag.concentrationPlot = "all";    % "all"/"last"/"none"
% file.Flag.potentialPlot = "all";
% file.Flag.currentPlot = "last";
% file.Flag.generationPlot = "none";
% file.Flag.characteristicCurvePlot = "Y";
% file.Flag.experimentalCurrentPotentialPlot = "Y";
% file.Flag.experimentalConcentrationPositiveIonPlot = "Y";
% file.Flag.experimentalPotentialPlot = "Y";
plotter(Res,file.Dati,file.Flag, file.Param);

