clear all;

f1Name = "alphaExpBeta7e5cells101.mat";
f2Name = "alphaExpBeta50kV_101cells.mat";
saveName = "alphaExpBeta50kV_101cells_full.mat";


f1 = load(fullfile(".\sim\", f1Name)).file;
f2 = load(fullfile(".\sim\", f2Name)).file;

%  fix for wrong adim
% f1xBarVec = [repmat(f1.ADati.Vbar, 1, f1.ADati.lr), repmat(f1.ADati.nbar, 1, f1.ADati.lr), repmat(f1.ADati.nbar, 1, f1.ADati.lr)]';
% f2xBarVec = [repmat(f2.ADati.Vbar, 1, f2.ADati.lr), repmat(f2.ADati.nbar, 1, f2.ADati.lr), repmat(f2.ADati.nbar, 1, f2.ADati.lr)]';
% f2.Res.Sol(:,1) = (f2.Res.Sol(:,1) ./ f2xBarVec) .* f1xBarVec;

file = f1;
file.Dati.T = f1.Dati.T + f2.Dati.T;
file.Dati.K = f1.Dati.K + f2.Dati.K;
file.Res.ASol = [f1.Res.ASol, f2.Res.ASol];
file.Res.Sol = [f1.Res.Sol, f2.Res.Sol];
file.Res.elapsedTime = f1.Res.elapsedTime + f2.Res.elapsedTime;
file.Res.kf = f1.Res.kf + f2.Res.kf;

save(fullfile(".\sim\", saveName), 'file');
fprintf("\nSaved Solution to %s.m \n",saveName);