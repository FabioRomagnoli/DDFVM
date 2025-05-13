clear all

% PLASMA CASES
% loadSol = "plasmaConstGen.mat";               % constant generation term in ion zone
% loadSol = "plasmaConstAlpha.mat";             % somehow ran to the end with const alpha but very wrong current
% loadSol = "alphaExpBeta3e5.full.mat";         % with beta found through const gen term comparison
% loadSol = "alphaExpBeta7e5.ful.mat";          % 1.61490e-04  first result that worked beta value from papers
% loadSol = "alphaExpBeta7e5newMu.full.mat";    % 1.608e-04, params (beta, Ei, mu) taken from papers,  51  cells, 2 sec
% loadSol = "alphaExpBeta7e5-100sec.mat";       % 1.60822e-04  100 sec 
% loadSol = "alphaExpBeta7e5cells101.mat";      % 1.36e-04 100 cells worst reslts


% DIODE CASES
% loadSol = "diode\diodeCase1.mat"
% loadSol = "diode\diodeCase2.mat"
% loadSol = "diode\diodeCase3.mat"
% loadSol = "diode\diodeSplit.mat"


% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s\n",loadSol);

% OUTPUT CONFIGURATON
printConfiguration(struct(),file.Param,struct(),file.Flag,struct(),struct());

% POST PROCESSING IF NECESSARY
file.Flag.saveSol = "no";       % IMPORTANT IT REMAINS "no"
file.Flag.compAlpha = false;
file.Flag.alpha = "const";  % "const"/"exp"
Res = postProcess(file.Dati, file.ADati, file.Res, file.Flag, file.Param);


% PLOTTING
file.Flag.concentrationPlot = "none";    % "all"/"last"/"none"
file.Flag.potentialPlot = "none";
file.Flag.currentPlot = "none";
file.Flag.generationPlot = "none";

plotter(file.Res,file.Dati,file.Flag);
