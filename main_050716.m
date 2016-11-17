addpath(genpath('/home/allison/Libraries/'))
% savepath allison@serenity:/home/allison/Libraries/pathdef.m

clear 
close all
% clc

%% - - - - - - - - - - LOAD & PROCESS RAW ENVIRONMENTS  - - - - - - - - - - %%
% % PROCESS RAW DATA
% raw2src_2_1_Tb  % temperature
% load('src_allTemporals.mat') 
% raw2src_2_1_X   % food availability
% load('srcChl_allTemporals.mat')

%% - - - - - - - - - - RUN OR LOAD BOOTSTRAP  - - - - - - - - - - %%
% run_Bootstrap_2_1   % temperature
load('bootData_w20Reps.mat')
% runChl_Bootstrap_2_1    % food
load('bootDataChl_w20Reps.mat')

%% - - - - - - - - - - SETUP DEB MODEL  - - - - - - - - - - %%
% categorize data
ftemporalNames = {'h';'d';'week';'mnth';'seasn';'yr'};
% ftemporalNames = {'h';'yr'};
temporalNames = {'h';'d';'7d';'30d';'120d';'365d'};
% temporalNames = {'h';'365d'};
spatialNames = {'BB';'BBSST';'LOBO';'off';'on';'SH';'SHSST'};
% spatialNames = {'off';'SH'};
foodNames = {'BBSSTc';'LOBOc';'SHSSTc'};
% foodNames = {'SHSSTc'};

% set sampling interval (hourly)
smpl_interval = 1/24;

% # replicates for ea DEB model run
nreps = 3;

% create time series time points starting with the first date/time in original raw
% data
t1 = datenum('25-Sep-2010 07:00:00'); ndtnum = t1:1/24:addtodate(t1,20,'year');
% t1 = datenum('25-Sep-2010 07:00:00'); ndtnum = t1:1/24:addtodate(t1,10,'year');

%% - - - - - - - - - - RUN DEB MODEL  - - - - - - - - - - %%
% following code replaces 'ndata2DEBout_wChldata.m'

% Model parameters (set outside DEB model function to save time):
% load and set structure with metaData, par data and auxiliary data (files saved from last estimation
% routine)
load('metaData.mat')
load('data.mat')
% load('auxData.mat')
estimData.metaData = metaData;  estimData.data = data;

% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);


for i = 1:length(temporalNames) % temporal
    tic
    for y = 1:length(spatialNames)
        src.(ftemporalNames{i})(y) = cellstr(strcat(spatialNames{y},'_',temporalNames{i}));
    end
    for z = 1:length(foodNames)
        srcChl.(ftemporalNames{i})(z) = cellstr(strcat(foodNames{z},'_',temporalNames{i}));
    end
%         srcNames = fieldnames(src.(temporalNames{i}));
        srcNames = src.(ftemporalNames{i});
        Xnames = srcChl.(ftemporalNames{i});
      
    for j = 1:length(srcNames)  % spatial
        
        for k = 1:length(foodNames)   % food
%                 Xnames = fieldnames(srcChl.(temporalNames{i}));   

            for L = 1:nreps;    % replicates
                Tbenv = ndata.(srcNames{j})(:,L);
                Xenv = ndataChl.(Xnames{k})(:,L);
                % set initial food density
                X_0 = Xenv(1);
                % set initial condition parameter stucture
                initpars = get_initpars(X_0, par);
                [t, f, Tb_out, Tbcorr_out, X_out, Lw, W_w, W_d, TRO, nrepyrs, time2mat, E, V, E_H, E_R, e] = run_DEB_042016(par, initpars, smpl_interval, Tbenv, Xenv, ndtnum, 0);
                [Edata(L).t, Edata(L).f, Edata(L).Tb_out, Edata(L).Tbcorr_out, Edata(L).X_out, ...
                    Edata(L).Lw, Edata(L).W_w, Edata(L).W_d, Edata(L).TRO, Edata(L).nrepyrs, ...
                    Edata(L).time2mat, Edata(L).E, Edata(L).V, Edata(L).E_H, Edata(L).E_R, ...
                    Edata(L).e] = run_DEB_042016(par, initpars, smpl_interval, Tbenv, Xenv, ndtnum, 0);

%                 DEBout.(srcNames{j}).(Xnames{k}).t{L} = t;
%                 DEBout.(srcNames{j}).(Xnames{k}).Lw{L} = Lw;
%                 DEBout.(srcNames{j}).(Xnames{k}).Ww{L} = W_w;
%                 DEBout.(srcNames{j}).(Xnames{k}).Wd{L} = W_d;
%                 DEBout.(srcNames{j}).(Xnames{k}).TRO{L} = TRO;
%                 DEBout.(srcNames{j}).(Xnames{k}).nrepyrs{L} = nrepyrs;
%                 DEBout.(srcNames{j}).(Xnames{k}).time2mat{L} = time2mat;
%                 savefig(strcat('DEBout',srcNames{j},Xnames{k},num2str(L),'.fig'));
            end
            DEBout.(ftemporalNames{i})(length(spatialNames)*(j-1)+k).name_Tb = spatialNames{j};
            DEBout.(ftemporalNames{i})(length(spatialNames)*(j-1)+k).name_X = foodNames{k};
            DEBout.(ftemporalNames{i})(length(spatialNames)*(j-1)+k).predict_data = Edata;
            toc
        end
        
%         Uncomment next 2 lines to save outputs. Comment to save time (but doesn't
%         store subsequent runs)
%           out = DEBout.(srcNames{j});
%           save(strcat(srcNames{j},'.mat'),'-struct','out')
    end
    
end

% save(strcat('DEBout',num2str(nreps),'yrs.mat') , 'DEBout')
save('DEBout.mat','DEBout')






