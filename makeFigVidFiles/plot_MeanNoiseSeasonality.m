% plot all mean*noise*sasonality combination
clear
close all

doSave = 0;

nyrs = 10;                          % # years
tInterval = 1;                      % sampling frequency per day eg tInterval = 1: 1 point per day, tInterval = 24: 24 points per day
rho = 1;                            % rho: cross-correlation of time series; default = 1
gamma = 0;                          % gamma: relationship between frequency and power: P(f)=1/f^gamma; -2 < gamma < 0: blue noise; 0 < gamma < 2: red noise
mu = 0;                             % mu: mean of each time series; default = 0
muPlus = [27 20 12 8];              % muPlus: add to mean
seasonal_amp = linspace(0,10,4);    % seasonal_amp: seasonal amplitudes
nreps = 10;                         % nreps: # replicates
nt = floor(365.25*nyrs*tInterval);  % number of time points
sd = [5 3 1 0];                             % sd: std. deviation within each time series (original simulations = 3)
            
% plot 0 seasonality all noise
figure
colormap Gray
hold on

for a = 1:length(muPlus)
    for b = 1:length(sd)
        env = phase_partnered_timeseries (nt, rho, gamma, sd(b), muPlus(a), doSave); 
        plot(1:nt, env, 'Color', [0 0 0]+(1/3)*(b-1));
        hold on
    end
end
