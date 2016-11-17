function getTPCfig(smpl_interval, doSave)
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
T = (-10:0.5:45) + 273.15;   
T_1 = par.T1;
Tpars = [par.T_A par.T_L par.T_H par.T_AL par.T_AH];
TPC = tempcorr(T,T_1,Tpars);

figure
plot(T-273.15,TPC,'LineWidth',3,'color','k')
xlabel('Body temperature')
ylabel('Performance')

xlim([-10 45])
ylim([0 1.2])

ax = gca; set(ax,'FontSize',12);

if doSave
    figname = sprintf('Mcal_PerformanceCurve.tiff');
    matfigname = sprintf('Mcal_PerformanceCurve.fig');
    
        saveas(gcf, figname)
        saveas(gcf, matfigname)
end
    