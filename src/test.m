clear all

% Residual Plot
loadSol = "diode\diodeCase2.mat";
sol1 = load(fullfile(".\sim\", loadSol)).file;

% 
% loadSol = "diode\diodeCase2.mat";
% sol2 =  load(fullfile(".\sim\", loadSol)).file;
% 


%% N concentration
figure();
% res = sol1.Res.Sol(Adaptive.Dati.nIdxs,end) - sol2.Res.Sol(sol1.Dati.nIdxs,end); 

% plot(res./sol1.Res.Sol(sol1.Dati.nIdxs,end));   % normalized residual of n  sol1 against sol2
% semilogy(sol1.Res.Sol(sol1.Dati.nIdxs,end))   % n solution at last time step of  sol1
% semilogy(sol2.Res.Sol(sol2.Dati.nIdxs,end))   % n sol2

figure();
hold on;
semilogy(sol1.Dati.r,sol1.Res.Sol(sol1.Dati.nIdxs,1),'DisplayName',"n")   % first concentration
semilogy(sol1.Dati.r,sol1.Res.Sol(sol1.Dati.pIdxs,1),'DisplayName',"p")   % first concentration
hold off;
set(gca, 'YScale', 'log') % Change y-axis to log scale
legend('Location', 'best'); 


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

