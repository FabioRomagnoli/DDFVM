clear all

% PLASMA CASES
loadSol = "alphaExpBeta7e5-100sec_50kV_full.mat";

% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s",loadSol);

% POST PROCESSING IF NECESSARY
file.Flag.saveSol = loadSol;     % !!Risk of overriding
file.Flag.computeCurrents = true;

postProcess(file.Dati, file.ADati, file.Res, file.Flag, file.Param);