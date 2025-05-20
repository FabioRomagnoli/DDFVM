clear all


% Residual Plot

loadSol = "diode\diodeConv51Adaptive.mat";
sol1 = load(fullfile(".\sim\", loadSol)).file;


loadSol = "diode\diodeConv51NotAdaptive.mat";
sol2 =  load(fullfile(".\sim\", loadSol)).file;



%% N concentration
figure();
res = sol1.Res.Sol(Adaptive.Dati.nIdxs,end) - sol2.Res.Sol(sol1.Dati.nIdxs,end); 

plot(res./sol1.Res.Sol(sol1.Dati.nIdxs,end));   % normalized residual of n  sol1 against sol2
% semilogy(sol1.Res.Sol(sol1.Dati.nIdxs,end))   % n solution at last time step of  sol1
% semilogy(sol2.Res.Sol(sol2.Dati.nIdxs,end))   % n sol2

% semilogy(sol1.Res.Sol(sol1.Dati.nIdxs,end))   % first concentration



%% P concentration
figure();
res = sol1.Res.Sol(Adaptive.Dati.pIdxs,end) - sol2.Res.Sol(sol1.Dati.pIdxs,end); 

plot(res./sol1.Res.Sol(sol1.Dati.pIdxs,end));   % normalized residual of p sol1 against sol2
% semilogy(sol1.Res.Sol(sol1.Dati.pIdxs,end))   % p solution at last time step of  sol1
% semilogy(sol2.Res.Sol(sol2.Dati.pIdxs,end))   % p sol2

% semilogy(sol1.Res.Sol(sol1.Dati.pIdxs,end))   % first concentration of  p




%% N concentration
figure();
res = sol1.Res.Sol(Adaptive.Dati.vIdxs,end) - sol2.Res.Sol(sol1.Dati.vIdxs,end); 

plot(res./sol1.Res.Sol(sol1.Dati.vIdxs,end));   % normalized residual of v  sol1 against sol2
% semilogy(sol1.Res.Sol(sol1.Dati.vIdxs,end))   % v solution at last time step of  sol1
% semilogy(sol2.Res.Sol(sol2.Dati.vIdxs,end))   % v sol2

% semilogy(sol1.Res.Sol(sol1.Dati.vIdxs,end))   % first concentration

