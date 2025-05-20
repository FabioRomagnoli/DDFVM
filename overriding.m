clear all

% PLASMA CASES
% loadSol = "plasmaConstGen.mat";               % constant generation term in ion zone
% loadSol = "plasmaConstAlpha.mat";             % very short sim time with const alpha but very wrong current
% loadSol = "alphaExpBeta3e5.full.mat";         % with beta found through const gen term comparison
% loadSol = "alphaExpBeta7e5.full.mat";         % 1.614e-04  first result that worked beta value from papers
loadSol = "alphaExpBeta7e5newMu.full.mat";    % 1.608e-04, params mu taken from papers 2 sec
% loadSol = "alphaExpBeta7e5-100sec.mat";       % 1.608e-04  100 sec 
% loadSol = "alphaExpBeta7e5cells101.mat";      % 1.360e-04 100 cells worst reslts

% LOAD SOLUTION AND DEF PARAMETERS
load(fullfile(".\sim\", loadSol));
fprintf("\nLoaded Solution %s",loadSol);

% POST PROCESSING IF NECESSARY
file.Flag.saveSol = loadSol;     % !!Risk of overriding
file.Flag.computeCurrents = true;

postProcess(file.Dati, file.ADati, file.Res, file.Flag, file.Param);