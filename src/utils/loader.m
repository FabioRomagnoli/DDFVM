clear all


% PLASMA CASES
% loadSol = "plasmaConstGen.mat"; % constant generation term in ion zone
% loadSol = "plasmaConstAlpha.mat"; % somehow ran to the end with const alpha but very wrong current
% loadSol = "alphaExpBeta7e5newMu.full.mat"; % 1.608e-04, params (beta, Ei, mu) taken from papers,  101  cells, 0.01 sec
% loadSol = "alphaExpBeta7e5cells101.mat"; % 1.36e-04
% loadSol = "alphaExpBeta7e5-100sec.mat";

loadSol = "alphaExpBeta3e5.full.mat";


% DIODE CASES
% loadSol = "dioceCase1.mat"
% loadSol = "dioceCase2.mat"
% loadSol = "diodeCase3.mat"
% loadSol = "diodeSplit.mat"


% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s\n",loadSol);
if strcmp(file.Flag.model,"plasma"); [Param, ~, ~] = initPlasma(struct(), file.Flag); end
if strcmp(file.Flag.model,"diode"); [Param, ~, ~] = initDiode(struct(), file.Flag); end

% OUTPUT CONFIGURATON
printConfiguration(file.Dati,Param,struct(),file.Flag,struct(),struct());

% POST PROCESSING IF NECESSARY
file.Flag.saveSol = "no";       % IMPORTANT IT REMAINS "no"
file.Flag.compAlpha = false;
file.Flag.alpha = "const";  % "const"/"exp"
Res = postProcess(file.Dati, file.ADati, file.Res, file.Flag, Param);


% PLOTTING
file.Flag.concentrationPlot = "last";    % "all"/"last"/"none"
file.Flag.potentialPlot = "last";
file.Flag.currentPlot = "last";
file.Flag.generationPlot = "none";

plotter(file.Res,file.Dati,file.Flag);
