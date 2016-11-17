% % % % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Part 1: create time series to use in DEB simulations
% Options:    run file: run Run phase_partnered_timeseries.m
%             # years
%             sampling interval: hourly, daily, etc
%             rho: (cross-correlation) unnecessary if just interested in making single 'base' time series. 
%             gamma: relationship between frequency and power: P(f)=1/f^gamma
%             sd: standard deviation
%             mu: mean; set to 0 b/c ultimately creating time series with different means. Faster to do it later.
%             muPlus: values to add to 0 mean
%             seasonal_amp: vector of amplitudes for seasonal components
%             doSave: save outputs to .mat files
%             doPlot: plot outputs
%             nreps: number of replicate time series to be created
% 
% Part 2: plot body temperature and relative performance PDFs
% Options:    doPlotTbw: 0 or 1
% 
% Part 3: run DEB model
% Options:    runDEB: 1 or 0
%             smpl_intervl: sampling interval; hourly, daily, etc. Must correspond
%             with tInterval (Part 1)
%             foodDensity: Set food densities
%
% Part 4: run analyze_correlated_environments
% Options:    See file
% % % % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clear
close all
clc

% runPhasePartnered: get phase partnered times series? yes = 1, no = 0
runPhasePartnered = 0;
    % doSave: save results? 1=yes, 0=no
    doSave = 0;
    % doPlot: plot results? 1=yes, 0=no
    doPlot = 0;
    
doPlotTbw = 0;  % plot Tb and relative performance PDF? yes = 1, no = 0

% runDEB: run DEB model sims? 1=yes, 0=no
runDEB = 0;

% runACE: run analyze_correlated_environments? yes = 1, no = 0
runACE = 1;

%% create base time series
% SET: nyrs, sd, muPlus, seasonal_amp, nreps, nt

nyrs = 10;                          % # years
tInterval = 1;                      % sampling frequency per day eg tInterval = 1: 1 point per day, tInterval = 24: 24 points per day
rho = 1;                            % rho: cross-correlation of time series; default = 1
gamma = 0;                          % gamma: relationship between frequency and power: P(f)=1/f^gamma; -2 < gamma < 0: blue noise; 0 < gamma < 2: red noise
sd = 10;                             % sd: std. deviation within each time series (original simulations = 3)
mu = 0;                             % mu: mean of each time series; default = 0
muPlus = [8 12 20 27];              % muPlus: add to mean
seasonal_amp = linspace(0,10,4);    % seasonal_amp: seasonal amplitudes
nreps = 10;                         % nreps: # replicates
nt = floor(365.25*nyrs*tInterval);  % number of time points

% create datenum vector (for plotting)
t1 = datenum('25-Sep-2010 07:00:00');   % day 1 of real temp data
dtnum = t1:1/tInterval:addtodate(t1,nyrs,'year'); dtnum = dtnum';

if(runPhasePartnered)    
    for k = 1:nreps
        env = phase_partnered_timeseries (nt, rho, gamma, sd, mu, doSave);
        env_0(:,k) = env;
    end
    env = env_0;
    if (doSave)
        fname = sprintf('environment-env_plus_0.000-env_samp_0.000.mat');
        save(fname, 'env');
    end
    if (doPlot)
        pname = sprintf('environment- env plus 0.000; env samp 0.000.mat');
        figure
        plot_ts(dtnum(1:size(env,1)), env, 'Body Temperature', pname)
    end
    
%% ----- 1.b adjust the mean
    if doPlot
        figure
    end
    load(fname)
    for k = 1:length(muPlus)
        envPlus = env + muPlus(k);
        if (doSave)
            fname = sprintf('environment-env_plus_%1.3f-env_samp_0.000.mat',muPlus(k));
            save(fname, 'envPlus');
        end
        if (doPlot)
            pname = sprintf('environment-env plus %1.3f, %1.3f, %1.3f, %1.3f; env samp 0.000',muPlus(1:4));
            plot_ts(dtnum(1:size(env,1)), envPlus, 'Body Temperature', pname)
            ylim([-2 42])
            hold on
        end
    end
    
%% ----- 1.c add seasonality
    time = 1:tInterval:nt;
    freq=1;
    sampling_freq=1/(nt/nyrs);
    
    for l = 1:length(seasonal_amp)
        s(:,l) = seasonal_amp(l) * cos(2*pi*freq*time*sampling_freq);
    end
    
    if doPlot
        figure
        hold on
    end
    for k = 1:length(muPlus)
        fname = sprintf('environment-env_plus_%1.3f-env_samp_0.000.mat', muPlus(k));
        load(fname)
        
        for l = 1:length(seasonal_amp)
            envSeason = bsxfun(@plus,envPlus,s(:,l));
            
            if (doSave)
                fname = sprintf('environment-env_plus_%1.3f-env_samp_%1.3f.mat', muPlus(k), seasonal_amp(l));
                save(fname, 'envSeason');
            end
            if (doPlot)
                pname = sprintf('environment-env plus %1.3f-env samp %1.3f', muPlus(k), seasonal_amp(l));
                subplot(length(seasonal_amp) , 1 , l)
                plot_ts(dtnum(1:size(env,1)), envSeason, 'Body Temperature', pname)
                ylim([-15 50])
                hold on
            end
        end
    end
end

%% % % % % ----- PART 2 ----- % % % % %
% Plot body temp and rel performance PDF
if doPlotTbw
    % extract data by filename organization
    Xlabel = 'Body temperature, °C';
    Ylabel = 'Probability density';
    
    colorOrder = [0 0.447 0.741;
        0.466 0.674 0.188;
        0.929 0.694 0.125;
        0.850 0.325 0.098];
    
    for x = 1:length(seasonal_amp)
        figure
        hold on
        for y = 1:length(muPlus)
            fname = sprintf('environment-env_plus_%1.3f-env_samp_%1.3f_10sdNoise.mat', muPlus(y), seasonal_amp(x));
            pname = sprintf('env plus %1.3f-env samp %1.3f', muPlus(y), seasonal_amp(x));
            load(fname)
            
            fdist = envSeason(:,1);
            [f,xi] = ksdensity(fdist,'Function','pdf');
            
            subplot(2,1,1)
            plot_BGYRhist(xi, f, y, Xlabel, Ylabel, pname);
            hold on
            
            subplot(2,1,2)
            TPC = getTPC(1,fdist);
            [f,xi] = ksdensity(TPC(TPC<= 1 & TPC >=0),'Function','pdf','Bandwidth',0.03);
            
            plot_BGYRhist(xi, f, y, 'Relative Performance', Ylabel, '');
            xlim([0 1])
            hold on
        end
    end
end

%% % % % % ----- PART 3 ----- % % % % %
    % set sampling interval: 1/24=hourly, 1=daily
    smpl_interval = 1;
    
    % set food values (here constant)
    % half saturation coefficient K = 6.12
    
    % based on percentiles
    %     chl_low = 3.55; % 33rd percentile of bootstrapped (n=20) GlobColour SH daiy, f = 0.367
    %     chl_med = 8.32; % 66th percentile, f = 0.576
    %     chl_hi = 44.58; % 99th percentile, f = 0.879
    % based on relative to f
    % X = getX(f,K)
    chl_low = 2.62; % f = 0.3
    %     chl_low = 4.08; % f = 0.4
    chl_med = 6.12; % f = 0.5
    chl_hi = 55.08; % f = 0.9
    
    % pack food densities
    foodDensity = [chl_low chl_med chl_hi];
% Run DEB model
if (runDEB)
    % Model parameters (set outside DEB model function to save time):
    % load and set structure with metaData, par data and auxiliary data (files saved from last estimation routine)
    load('metaData.mat')
    load('data.mat')
    % load('auxData.mat')
    estimData.metaData = metaData;  estimData.data = data;
    
    % set par structure
    [par, txtPar] = set_pars(smpl_interval, estimData);
    
    for m = 1:length(foodDensity) % # of food densities
        Xenv = repmat(foodDensity(m),nt,1);
        tic
        % set initial condition parameter stucture (food dependent)
        initpars = get_initpars(foodDensity(m), par);
        
        for k = 1:length(muPlus)
            fname = sprintf('environment-env_plus_%1.3f-env_samp_0.000_10sdNoise.mat', muPlus(k));
            load(fname)
            for l = 1:length(seasonal_amp)
                fname = sprintf('environment-env_plus_%1.3f-env_samp_%1.3f_10sdNoise.mat', muPlus(k), seasonal_amp(l));
                load(fname)
                for n = 1:nreps
                    Tbenv = envSeason(:,n);
                    out = run_DEB_071116(par, initpars, smpl_interval, Tbenv, Xenv, dtnum, 0);
                    DEB_out(n) = out;
                end
                fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_10sdNoise.mat', muPlus(k), seasonal_amp(l), foodDensity(m));
                save(fname,'DEB_out')
            end
        end
    end
end

%% run analyze_correlated_environment. Contains ScatteredInterpolation, compare_noise, etc.
% names = analyze_correlated_environment(gamma, muPlus, seasonal_amp);
if runACE
analyze_correlated_environment
end

