% plot DEB outputs across seasonal components and food levels
% ea plot has 4 mean temperatures

plotsx = 0;

sd = [0 1 3 5 10];


% plotDEBout
plotDEBout = 0;
if plotDEBout
    for z = 1:length(foodDensity)
        for y = 1:length(seasonal_amp)
            figure
            hold on
            
            for x = 1:length(muPlus)
                fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', muPlus(x), seasonal_amp(y), foodDensity(z));
                load(fname)
                data = DEB_out(1); mainTitle = '';
                Title = strcat('+', int2str(muPlus), ' samp=', int2str(seasonal_amp(y)), ' X=', int2str(foodDensity(z)));
                get_xPlots
                %                 figure
                %                 histogram(DEB_out(1).Tb_out)
                %                 hold on
            end
            sp1 = mtit(Title, 'fontsize', 12, 'xoff',-.1,'yoff',.025);
            for i = 1:length(muPlus)
                Legend{i} = strcat('mu = ', int2str(muPlus(i)));
            end
            legend(Legend,'Location','best');
        end
    end
end

% plot ingestion rate, assimilation rate, growth
plotpXpApG = 0;

% load parameters
load('metaData.mat')
load('data.mat')
estimData.metaData = metaData;  estimData.data = data;
% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);
cPar = parscomp_st(par);
vars_pull(cPar);

% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);
if plotpXpApG
    for z = 1:length(foodDensity)
        
        for y = 1:length(seasonal_amp)
            figure
            hold on
            
            for x = 1:length(muPlus)
                fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', muPlus(x), seasonal_amp(y), foodDensity(z));
                load(fname)
                data = DEB_out(1); mainTitle = '';
                
                f = data.f;
                E = data.E;
                V = data.V;
                T = data.Tb_out;
                X = data.X_out;
                Lw = data.Lw;
                
                initpars = get_initpars(X(1), par);
                M = shapecorr(Lw, initpars.L_b, initpars.L_j); % shape correction factor
                
                
                % TP parameters
                T_1 = par.T1;
                Tpars = [par.T_A par.T_L par.T_H par.T_AL par.T_AH];
                c_T = tempcorr(T+273.15,T_1,Tpars);
                
                % correct for temperature and shape
                p_AmT = c_T .* cPar.p_Am .* M; % specific maximum assimilation rate
                p_XmT = p_AmT ./ par.kap_X;  % max specific ingestion rate
                p_MT = c_T .* par.p_M;
                v_T = c_T .* par.v .* M;          % energy conductance
                
                % get rates
                pX = f .* p_XmT .* V.^(2/3);  % ingestion rate
                pA = par.kap_X .* pX;    % assmilation rate
                pM = p_MT .* V;             % somatic maintenance flux
                pS = pM;	% volume linked somatic maintenance, pS = pM + pT pT = 0 for bvalves)
                %                 pT = p_TT .* V.^(2/3);   % surface area-linked somatic maintenance; maintenance of osmotic regulatory boundaries
                pC = (E./V) .* (par.E_G .* v_T .* V.^(2/3) + pS ) ./ (par.kap .* E./V + par.E_G ); %eq. 2.12 p.37 Kooijman 2010
                pG = (par.kap .* pC - pS)/ par.E_G;        % growth (in structural volume)
                G = pG.^(1/3); % conert to structural length
                Gw = G ./ par.del_M;
                
                for i = 1:length(data.W_w)-1
                    GR(i) = (real(data.W_w(i+1)) - real(data.W_w(i)))./real(data.W_w(i));
                end
                
                
                Title = strcat('+', int2str(muPlus), ' samp=', int2str(seasonal_amp(y)), ' X=', int2str(foodDensity(z)));
                
                subplot(2, 1, 1)
                plot(data.t, pX, 'LineWidth', 1)
                hold on
                xlabel('Date')
                ylabel('-')
                title('ingestion rate')
                
                %                 subplot(4, y+y-1, 2)
                %                 plot(data.t, pA, 'LineWidth', 1)
                %                 hold on
                %                 xlabel('Date')
                %                 ylabel('-')
                %                 title('Assimilation rate')
                
                %                 subplot(4, y+y-1, 3)
                %                 plot(data.t, pM, 'LineWidth', 1)
                %                 hold on
                %                 xlabel('Date')
                %                 ylabel('-')
                %                 title('Somaic maintenance flux')
                
                subplot(2, 1, 2)
                plot(data.t(1:end-1), GR(1:length(data.t)-1), 'LineWidth', 1)
                hold on
                xlabel('Date')
                ylabel('-')
                title('Growth')
                
            end
            sp1 = mtit(Title, 'fontsize', 12, 'xoff',-.1,'yoff',.025);
            for i = 1:length(muPlus)
                Legend{i} = strcat('mu = ', int2str(muPlus(i)));
            end
            legend(Legend,'Location','best');
        end
    end
end

n_t = 365.25*5;
nnt = 1:n_t;

[map,~,~] = brewermap(11,'*Spectral');
% define colors for current and total TP
muPlusof1_TP = map(1,:);    muPlusof1_pt = map(2,:);
muPlusof2_TP = map(4,:);    muPlusof2_pt = map(5,:);
muPlusof3_TP = map(8,:);    muPlusof3_pt = map(7,:);
muPlusof4_TP = map(11,:);    muPlusof4_pt = map(10,:);

nmap = [muPlusof1_TP muPlusof1_pt; ...
    muPlusof2_TP muPlusof2_pt; ...
    muPlusof3_TP muPlusof3_pt; ...
    muPlusof4_TP muPlusof4_pt];
i_F = zeros(10,4);


if plotsx
    for z = 1:length(foodDensity)
        for y = 1:length(seasonal_amp)
            for w = 1:length(sd)
                for x = 1:length(muPlus)
                    fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_%dsdNoise.mat', muPlus(x), seasonal_amp(y), foodDensity(z), sd(w));
                    load(fname)
                    data = DEB_out(1); mainTitle = '';
                    Tb = interp1(data.t, data.Tb_out, 1:nt);
                    T_b(:,x) = Tb;
                    
                    %% comment
                    %             makeTbVid(seasonal_amp, y, foodDensity, z, nnt, Tb, muPlus, x, nmap(x,:))
                    %             makeTPCVid(seasonal_amp, y, foodDensity, z, nnt, Tb, muPlus, x, nmap(x,:))
                    %             figname = sprintf('TPC-env_plus_%1.3f-env_samp_%1.3f', muPlus(x), seasonal_amp(y));
                    %             makeTPCfig(seasonal_amp, y, Tb, muPlus, x, nmap(x,:));
                    %             saveas(gcf, figname)
                    %
                    %             plotTb_on_TPC(seasonal_amp, y, Tb, muPlus, x)
                    %             figname = strcat('TbonTPC_',int2str(muPlus(x)), '_', int2str(seasonal_amp(y)));
                    %             figname = char(figname);
                    %             saveas(gcf, figname, 'svg')
                    
                    %% use
                    % get DEB data
                    t = data.t;
                    f = data.f;
                    E = data.E;
                    V = data.V;
                    T = data.Tb_out;
                    X = data.X_out;
                    Lw = data.Lw;
                    Lw2 = interp1(data.t, Lw, 1:nt);
                    L_w(:,x) = Lw2;
                    Ww = data.W_w;
                    Ww = interp1(data.t, Ww, 1:nt);
                    W_w(:,x) = Ww;
                    
                    TRO = data.TRO;
                    TR_O(y,x) = real(TRO);
                    nrep = data.nrepyrs;
                    avgRepro(y,x) = TRO/nrep;
                    
                    E_Hp = 40.78; kap_R = 0.95;   % FROM SET_PARS_NEW.M
                    % Fraw = zeros(size(data.rawE_R,1),size(data.rawE_R,2));
                    E_0 = E(1);
                    E_R = data.E_R;
                    E_H = data.E_H;
                    i_sp = find(and((E_R == 0),(E_H>=E_Hp)));
                    i_sp = i_sp - 1; % look at preceeding line for E_R value before spawning
                    F = kap_R .* E_R ./ E_0;  % fecundity = egg number
                    
                    i_spF = F(i_sp);
                    i1 = 10-length(i_spF)+1;
                    i_F(i1:end,x) = i_spF;
                    
                    initpars = get_initpars(X(1), par);
                    M = shapecorr(Lw, initpars.L_b, initpars.L_j); % shape correction factor
                    
                    % TP parameters
                    T_1 = par.T1;
                    Tpars = [par.T_A par.T_L par.T_H par.T_AL par.T_AH];
                    c_T = tempcorr(T+273.15,T_1,Tpars);
                    
                    % correct for temperature and shape
                    p_AmT = c_T .* cPar.p_Am .* M; % specific maximum assimilation rate
                    p_XmT = p_AmT ./ par.kap_X;  % max specific ingestion rate
                    
                    % get rates
                    pX = f .* p_XmT .* V.^(2/3);  % ingestion rate
                    pX = interp1(data.t, pX, 1:nt);
                    p_X(:,x) = pX;
                end
                %% make vid figs
                %         makepXvid(seasonal_amp, y, foodDensity, z, nnt, p_X, muPlus, nmap)
                %         makepXfig(seasonal_amp, y, foodDensity, z, nnt(365:365*2), p_X(365:365*2,:), muPlus, nmap)
                %         pXfigname = sprintf('pX-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(z));
                %         saveas(gcf, pXfigname)
                %
                %         makeTbfig(seasonal_amp, y, foodDensity, z, 1:nt, T_b, muPlus, nmap)
                %         Tbfigname = sprintf('Tb-all_muPlus-samp_%1.3f', seasonal_amp(y));
                % %         saveas(gcf, Tbfigname)
                %         saveas(gcf,Tbfigname, 'svg')
                
                %         makeWwfig(seasonal_amp, y, foodDensity, z, 1:nt, W_w, muPlus, nmap)
                %         Wwfigname = sprintf('Ww-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(z));
                %         saveas(gcf, Wwfigname)
                
                %         makeWwfig(seasonal_amp, y, foodDensity, z, 1:nt, L_w, muPlus, nmap)
                %         Lwfigname = sprintf('Lw-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(z));
                %         saveas(gcf, Lwfigname)
                
                %         makei_Ffig(seasonal_amp, y, foodDensity, z, nt, i_F, muPlus, nmap)
                %         Ffigname = sprintf('F-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(z));
                %         saveas(gcf, Ffigname)
                
                %% keep
                Wwmat(:,z,y) = real(max(W_w));
                Lwmat(:,z,y) = real(max(L_w));
                
                avgRepro = real(avgRepro);
                %     plot(y,avgRepro(y,:),'o')
            end
            
            %% comment
            % makeFfood_season_tempfig(avgRepro, seasonal_amp, foodDensity, muPlus, nmap)
            % s_Ffigname = sprintf('sxF-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(z));
            % saveas(gcf, s_Ffigname)
            
            % makeFfood_season_tempfig(TR_O, seasonal_amp, foodDensity, muPlus, nmap)
            % TROfigname = sprintf('TRO-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(z));
            % saveas(gcf, TROfigname)
        end
    end
    % path(p)
    Wwmat = real(Wwmat);
end



%%%%%
Figs = [1 0];   % [plotFigs saveFigs]
% ScatteredInterpolation3D(muPlus, seasonal_amp, foodDensity, sd, smpl_interval, Figs)
%%%%%
% MEScalcs(seasonal_amp, muPlus, foodDensity)
%%%%%
compare_noiseNonoise(seasonal_amp, muPlus, foodDensity, smpl_interval, sd)
%%%%%




% Lwmat = real(Lwmat);
% makefood_season_tempfig(Wwmat, seasonal_amp, foodDensity,muPlus,nmap)
% s_wwfigname = sprintf('sxWw-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(3));
%         saveas(gcf, s_wwfigname)
        
Fmat = real(Wwmat);
makefood_season_tempfig(Wwmat, seasonal_amp, foodDensity,muPlus,nmap)
s_wwfigname = sprintf('sxWw-env_plus_%1.3f%1.3f%1.3f%1.3f-env_samp_%1.3f-env_X_%1.3f.tiff', muPlus(1:4), seasonal_amp(y), foodDensity(3));
        saveas(gcf, s_wwfigname)

% makeffig(foodDensity)
% Xfigname = sprintf('Food Density-Functional Response.tiff');
% saveas(gcf, Xfigname)

disp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTb_on_TPC(seasonal_amp, y, T_b, muPlus, x)
TPC = getTPC(1, -5:40);

figure
plot(-5:40, TPC, 'k')
hold on

minTb = min(T_b); maxTb = max(T_b);
Tb = minTb:1:maxTb;
TbTPC = getTPC(1,Tb);
plot(Tb, TbTPC, 'r', 'LineWidth', 3)

title(strcat('mu Tb = ', int2str(muPlus(x)), ' samp=', int2str(seasonal_amp(y))));
end
