% Create two time series with specific autocorrelation and rho cross-correlation 
% using the algorithm described by Vasseur (2007)
% 
% [env] = phase_partnered_timeseries (n, rho, gamma, sd, mean, doSave)
%
%    n       = number of time steps
%    rho     = cross-correlation of time series
%    gamma   = relationship between frequency and power: P(f)=1/f^gamma
%              (-2 < gamma < 0: blue noise; 0 < gamma < 2: red noise
%    sd      = std. deviation within each time series
%    mean    = mean of each time series
%    doSave  = save results? 1=yes, 0=no
%    env     = n x 2 matrix containing the two time series

% Author: Tarik C. Gouhier
% Last modified: December 15, 2008

function [env] = phase_partnered_timeseries (n, rho, gamma, sd, mean, doSave)
if rho == 1
    env = ones(n,1)*NaN;
else
    env=ones(n,2)*NaN;
end
minVal=0;
maxVal=2*pi;
fs=[1:n/2]';
epsilon=rand(n/2,1);
epsilon (epsilon > 0.5)=1;
epsilon (epsilon < 0.5)=-1;

delta=acos(rho).*epsilon;
phi1=minVal+(maxVal-minVal)*rand(n/2,1);
phi2=phi1+delta;

if rho == 1
    for t=1:n
        env(t,1)=sum((1./(fs.^(gamma/2))).*(sin(2*pi.*fs*t/n + phi1)));
    end
else
    
    for t=1:n
        env(t,1)=sum((1./(fs.^(gamma/2))).*(sin(2*pi.*fs*t/n + phi1)));
        env(t,2)=sum((1./(fs.^(gamma/2))).*(sin(2*pi.*fs*t/n + phi2)));
    end
end

env=env./(repmat(std(env), [n, 1]));
env=mean+env.*sd;

% if (doSave)
%     fname=sprintf('correlated_environment_env_std_%1.3f-env_corr_%1.3f-gamma_%1.3f.txt', ...
%         sd, rho, gamma); 
%    save(fname, 'env', '-ascii'); 
% end
