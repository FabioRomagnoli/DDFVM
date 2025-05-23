% clear all

% PLASMA CASES
% loadSol = "plasmaConstGen.mat";               % constant generation term in ion zone
% loadSol = "plasmaConstAlpha.mat";             % very short sim time with const alpha but very wrong current
% loadSol = "alphaExpBeta7e5cells101.mat";      % 1.360e-04 100 cells worst reslts

% loadSol = "alphaExpBeta3e5.full.mat";         % with beta found through const gen term comparison
% loadSol = "alphaExpBeta7e5.full.mat";         % 1.614e-04  first result that worked beta value from papers
% loadSol = "alphaExpBeta7e5newMu.full.mat";    % 1.608e-04, params mu taken from papers 2 sec
% loadSol = "alphaExpBeta7e5-100sec.mat";       % 1.608e-04  100 sec 
% loadSol = "plasma\alphaExpBeta3e5.full.mat";      
% loadSol = "alphaExpBeta7e5-100sec_50kV_full.mat";
% loadSol = "alphaExpBeta50kV_81cells.mat";

% DIODE CASES
% loadSol = "diode\diodeCase1.mat";         % Given in Exam 
% loadSol = "diode\diodeCase2.mat";         % Given in Exam, assymmetricl
% loadSol = "diode\diodeCase3.mat";
% loadSol = "diode\diodeCase4.mat";
% loadSol = "diode\diodeCase5.mat";

% loadSol = "diode\diodeCase3inTime.mat";
% loadSol = "diode\diodeSplit.mat";
% loadSol = "diode\diodeReverse.mat";
loadSol = "diode\diodeSweep.mat";
% loadSol = "checkpoint.mat";

% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s\n",loadSol);

% OUTPUT CONFIGURATON
printConfiguration(struct(),file.Param,struct(),file.Flag,struct(),struct());


% % POST PROCESSING IF NECESSARY
file.Flag.saveSol = "no";     % !!Risk of overriding   /   loadSol
file.Flag.computeCurrents = true;
Res = postProcess(file.Dati, file.ADati, file.Res, file.Flag, file.Param);


% PLOTTING
% file.Flag.concentrationPlot = "reduced";    % "all"/"last"/"none"
% file.Flag.potentialPlot = "reduced";
% file.Flag.currentPlot = "reduced";
% file.Flag.generationPlot = "none";
file.Flag.characteristicCurvePlot = "Y";
% file.Flag.experimentalCurrentPotentialPlot = "Y";
% file.Flag.experimentalConcentrationPositiveIonPlot = "Y";
% file.Flag.experimentalPotentialPlot = "Y";
plotter(Res,file.Dati,file.Flag, file.Param);

