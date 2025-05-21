clear all;
dataPlourabu =  csvread("dati_esperimento_zheng\xx_PlourabuèData.csv");
% dataPlourabu(all(dataPlourabu == 0, 2), :) = [];  % removes rows where both elements are zero



dataNodes =  csvread("dati_esperimento_zheng\nodes.csv");
dataTime =  csvread("dati_esperimento_zheng\time_instants.csv");

dataPositiveIon =  csvread("dati_esperimento_zheng\positive.csv");
dataPotential =  csvread("dati_esperimento_zheng\potential.csv");


plot(dataNodes, dataPositiveIon(end,:));

% plot(dataNodes, dataPotential(end,:));

