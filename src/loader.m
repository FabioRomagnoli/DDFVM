% clear all

% PLASMA CASES
% loadSol = "plasmaConstGen.mat";               % constant generation term in ion zone
% loadSol = "plasmaConstAlpha.mat";             % very short sim time with const alpha but very wrong current
% loadSol = "alphaExpBeta7e5cells101.mat";      % 1.360e-04 100 cells worst reslts


% loadSol = "alphaExpBeta3e5.full.mat";         % with beta found through const gen term comparison
% loadSol = "alphaExpBeta7e5.full.mat";         % 1.614e-04  first result that worked beta value from papers
% loadSol = "alphaExpBeta7e5newMu.full.mat";    % 1.608e-04, params mu taken from papers 2 sec
% loadSol = "alphaExpBeta7e5-100sec.mat";       % 1.608e-04  100 sec 
% loadSol = "alphaExpBeta7e5-100sec_50kV.mat";      
% loadSol = "alphaExpBeta7e5-100sec_50kV_full.mat";
% loadSol = "alphaExpBeta50kV_81cells.mat";

% DIODE CASES
% loadSol = "diode\diodeCase1.mat";
loadSol = "diode\diodeCase2.mat";
% loadSol = "diode\diodeCase3.mat";
% loadSol = "diode\diodeCase4.mat";
% loadSol = "diode\diodeCase5.mat";

% loadSol = "diode\diodeCase3inTime.mat";
% loadSol = "diode\diodeSplit.mat";
% loadSol = "diode\diodeReverse.mat";

% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s\n",loadSol);

% OUTPUT CONFIGURATON
printConfiguration(struct(),file.Param,struct(),file.Flag,struct(),struct());

% PLOTTING
file.Flag.concentrationPlot = "none";    % "all"/"last"/"none"
file.Flag.potentialPlot = "none";
file.Flag.currentPlot = "none";
% file.Flag.generationPlot = "none";
% file.Flag.characteristicCurvePlot = "Y";
% file.Flag.experimentalCurrentPotentialPlot = "Y";
% file.Flag.experimentalConcentrationPositiveIonPlot = "Y";
% file.Flag.experimentalPotentialPlot = "Y";
plotter(file.Res,file.Dati,file.Flag);

