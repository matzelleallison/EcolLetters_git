
function TPC = getTPC(smpl_interval,T)
% Model parameters (set outside DEB model function to save time):
% load and set structure with metaData, par data and auxiliary data (files saved from last estimation
% routine)
load('metaData.mat')
load('data.mat')
% load('auxData.mat')
estimData.metaData = metaData;  estimData.data = data;

% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);

% get performane curve
T = T + 273.15;   
T_1 = par.T1;
Tpars = [par.T_A par.T_L par.T_H par.T_AL par.T_AH];
TPC = tempcorr(T,T_1,Tpars);


% figure
% plot(T-273.15,TPC,'LineWidth',3,'color','k')
% xlabel('Body temperature')
% ylabel('Performance')