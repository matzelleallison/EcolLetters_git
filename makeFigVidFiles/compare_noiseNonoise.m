function compare_noiseNonoise(seasonal_amp, muPlus, foodDensity, smpl_interval, sd)
% function extracts, plots, and saves plots of (i) Lw, Ww, Wd and (ii) pX, pM, pWw,
% pLw
% ea figure has 4 simulations for ea figure (0 noise, 1sd noise, 3sd noise, 5sd
% noise)
% 1 mean Tb, food level and seasonal amplitude per figure
% called by: analyze_correlated_environment.m
% calls" getRates.m

plotrates = 0;
plotTbTPCfreq = 1;
plotTbTPCspline = 0;
TPprobability = 0;
Tbprobability = 1;

for l = 1:length(foodDensity)
    %     fd = repmat(foodDensity(l), 10, 1);
    for k = 1:length(seasonal_amp)
        %         sa = repmat(seasonal_amp(k), 10, 1);
        for j = 1:length(muPlus)
            %             Tb = repmat(muPlus(j), 10, 1);
            
            % % % extract sims with 10sd noise
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_10sdNoise.mat', ...
                muPlus(j), seasonal_amp(k), foodDensity(l));
            load(fname)
            data_10sdNoise = DEB_out(1);
            
            [Lw_10sd, Ww_10sd, Wd_10sd, pX_10sd, pM_10sd, pLw_10sd, pWw_10sd, Tb_out_10sd] = getRates(data_10sdNoise, smpl_interval);
            
            % % % extract sims with 5sd noise
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_5sdNoise.mat', ...
                muPlus(j), seasonal_amp(k), foodDensity(l));
            load(fname)
            data_5sdNoise = DEB_out(1);
            
            [Lw_5sd, Ww_5sd, Wd_5sd, pX_5sd, pM_5sd, pLw_5sd, pWw_5sd, Tb_out_5sd] = getRates(data_5sdNoise, smpl_interval);
            
            % % % extract sims with 3sd noise
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_3sdNoise.mat', ...
                muPlus(j), seasonal_amp(k), foodDensity(l));
            load(fname)
            data_3sdNoise = DEB_out(1);
            
            [Lw_3sd, Ww_3sd, Wd_3sd, pX_3sd, pM_3sd, pLw_3sd, pWw_3sd, Tb_out_3sd] = getRates(data_3sdNoise, smpl_interval);
            
            % % % extract sims with 1sd noise
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_1sdNoise.mat', ...
                muPlus(j), seasonal_amp(k), foodDensity(l));
            load(fname)
            data_1sdNoise = DEB_out(1);
            
            [Lw_1sd, Ww_1sd, Wd_1sd, pX_1sd, pM_1sd, pLw_1sd, pWw_1sd, Tb_out_1sd] = getRates(data_1sdNoise, smpl_interval);
            
            % % % extract sims with no noise
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_0sdNoise.mat', ...
                muPlus(j), seasonal_amp(k), foodDensity(l));
            load(fname)
            data_noNoise = DEB_out(1);
            
            [Lw_noNoise, Ww_noNoise, Wd_noNoise, pX_noNoise, pM_noNoise, pLw_noNoise, pWw_noNoise, Tb_out_0sdNoise] = getRates(data_noNoise, smpl_interval);
            
            %% get mean temperature and mean performance and mean max Wd
            
            mu_0sd = nanmean(Tb_out_0sdNoise);
            mu_1sd = nanmean(Tb_out_1sd);
            mu_3sd = nanmean(Tb_out_3sd);
            mu_5sd = nanmean(Tb_out_5sd);
            
            meanTb_tab{j,k} = [mu_0sd; mu_1sd; mu_3sd; mu_5sd];
            
            tp_noN = nanmean(getTPC(1,Tb_out_0sdNoise));
            tp_1sd = nanmean(getTPC(1,Tb_out_1sd));
            tp_3sd = nanmean(getTPC(1,Tb_out_3sd));
            tp_5sd = nanmean(getTPC(1,Tb_out_5sd));
            
            meanTP_tab{j,k} = [tp_noN; tp_1sd; tp_3sd; tp_5sd];
            
            max_noN = nanmax(Wd_noNoise);
            max_1sd = nanmax(Wd_1sd);
            max_3sd = nanmax(Wd_3sd);
            max_5sd = nanmax(Wd_5sd);
            
            maxWd_tab{j,k} = [max_noN; max_1sd; max_3sd; max_5sd];
            
            
            if plotrates
                %% plot rates noise vs. no noise
                mainTitle = sprintf('Tb %1.3f; amp %1.3f; food %1.3f', muPlus(j), seasonal_amp(k), foodDensity(l));
                figname = sprintf('pX-pM-pLw-pWw_Tb_%1.3f-amp_%1.3f-X_%1.3f_Noise-vs-noNoise.tiff', muPlus(j), seasonal_amp(k), foodDensity(l));
                matfigname = sprintf('pX-pM-pLw-pWw_Tb_%1.3f-amp_%1.3f-X_%1.3f_Noise-vs-noNoise.fig', muPlus(j), seasonal_amp(k), foodDensity(l));
                nt = 1 : (365.25*10);
                
                figure
                
                subplot(2,2,1)
                plot(nt, pX_5sd)
                hold on
                plot(nt, pX_3sd)
                plot(nt, pX_1sd)
                plot(nt, pX_noNoise, 'LineWidth', 2)
                xlabel('day'); ylabel('ingestion rate');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                subplot(2,2,2)
                plot(nt, pM_5sd)
                hold on
                plot(nt, pM_3sd)
                plot(nt, pM_1sd)
                plot(nt, pM_noNoise, 'LineWidth', 2)
                xlabel('day'); ylabel('somatic maintenance');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                subplot(2,2,3)
                plot(nt(2:end), pLw_5sd)
                hold on
                plot(nt(2:end), pLw_3sd)
                plot(nt(2:end), pLw_1sd)
                plot(nt(2:end), pLw_noNoise, 'LineWidth', 2)
                xlabel('day'); ylabel('growth rate (dLw/dt)');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                subplot(2,2,4)
                plot(nt(2:end), pWw_5sd)
                hold on
                plot(nt(2:end), pWw_3sd)
                plot(nt(2:end), pWw_1sd)
                plot(nt(2:end), pWw_noNoise, 'LineWidth', 2)
                xlabel('day'); ylabel('growth rate (dWw/dt)');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                tits = mtit(mainTitle, 'fontsize', 10);
                %             saveas(gcf, figname)
                %             saveas(gcf,matfigname)
                
                % % % plot lengths and weights noise vs. no noise
                nmainTitle = sprintf('Tb %1.3f; amp %1.3f; food %1.3f', muPlus(j), seasonal_amp(k), foodDensity(l));
                nfigname = sprintf('Lw-Ww-Wd_Tb_%1.3f-amp_%1.3f-X_%1.3f_Noise-vs-noNoise.tiff', muPlus(j), seasonal_amp(k), foodDensity(l));
                nmatfigname = sprintf('Lw-Ww-Wd_Tb_%1.3f-amp_%1.3f-X_%1.3f_Noise-vs-noNoise.fig', muPlus(j), seasonal_amp(k), foodDensity(l));
                nt = 1 : (365.25*10);
                
                figure
                
                subplot(1,3,1)
                plot(nt, Lw_5sd)
                hold on
                plot(nt, Lw_3sd)
                plot(nt, Lw_1sd)
                plot(nt, Lw_noNoise, 'LineWidth', 1)
                xlabel('day'); ylabel('length');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                subplot(1,3,2)
                plot(nt, Ww_5sd)
                hold on
                plot(nt, Ww_3sd)
                plot(nt, Ww_1sd)
                plot(nt, Ww_noNoise, 'LineWidth', 1)
                xlabel('day'); ylabel('wet weight');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                subplot(1,3,3)
                plot(nt, Wd_5sd)
                hold on
                plot(nt, Wd_3sd)
                plot(nt, Wd_1sd)
                plot(nt, Wd_noNoise, 'LineWidth', 1)
                xlabel('day'); ylabel('dry weight');
                Legend = {'5 sd','3 sd', '1 sd', 'no noise'};
                legend(Legend, 'Location', 'NorthWest')
                
                tits = mtit(nmainTitle, 'fontsize', 10);
                %             saveas(gcf, nfigname)
                %             saveas(gcf,nmatfigname)
            end
            
            %% plot histogram and overlay PDF
            if plotTbTPCfreq
                %     Tbedge = [
                tpedge = [0:0.03:1];
                if Tbprobability
                    %% 0 sd
                    figure
                    histogram(Tb_out_0sdNoise, 'Normalization', 'pdf')
                    hold on
                    pd = fitdist(Tb_out_0sdNoise', 'Normal');
                    x = xlim;
                    x = x(1):0.1:x(2);
                    y = pdf(pd,x);
                    plot(x,y,'LineWidth',2)
                    xlabel('body temperature')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    
                    
                    mTitle = sprintf('Tb %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(1));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(1));
                    matfigname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(1));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 1 sd
                    figure
                    histogram(Tb_out_1sd,'Normalization', 'pdf')
                    hold on
                    pd = fitdist(Tb_out_1sd', 'Normal');
                    x = xlim;
                    x = x(1):0.1:x(2);
                    y = pdf(pd,x);
                    plot(x,y,'LineWidth',2)
                    xlabel('body temperature')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    
                    
                    mTitle = sprintf('Tb %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(2));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(2));
                    matfigname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(2));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 3 sd
                    figure
                    histogram(Tb_out_3sd, 'Normalization', 'pdf')
                    hold on
                    pd = fitdist(Tb_out_3sd', 'Normal');
                    x = xlim;
                    x = x(1):0.1:x(2);
                    y = pdf(pd,x);
                    plot(x,y,'LineWidth',2)
                    xlabel('body temperature')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    
                    mTitle = sprintf('Tb %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(3));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(3));
                    matfigname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(3));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 5 sd
                    figure
                    histogram(Tb_out_5sd, 'Normalization', 'pdf')
                    hold on
                    pd = fitdist(Tb_out_5sd', 'Normal');
                    x = xlim;
                    x = x(1):0.1:x(2);
                    y = pdf(pd,x);
                    plot(x,y,'LineWidth',2)
                    xlabel('body temperature')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    
                    mTitle = sprintf('Tb %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(4));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(4));
                    matfigname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(4));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 10 sd
                    figure
                    histogram(Tb_out_10sd, 'Normalization', 'pdf')
                    hold on
                    pd = fitdist(Tb_out_10sd', 'Normal');
                    x = xlim;
                    x = x(1):0.1:x(2);
                    y = pdf(pd,x);
                    plot(x,y,'LineWidth',2)
                    xlabel('body temperature')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    
                    
                    mTitle = sprintf('Tb %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(5));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(5));
                    matfigname = sprintf('TbPDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(5));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                end
                %% TP probability
                %% 0sd
                if TPprobability
                    figure
                    tp = getTPC(1,Tb_out_0sdNoise);
                    histogram(tp, tpedge, 'Normalization',  'probability');
                    ylim([0 0.4])
                    hold on
                    %             pd = fitdist(tp', 'Normal');
                    %             x = xlim;
                    %             x = x(1):0.1:x(2);
                    %             y = pdf(pd,x);
                    %             plot(x,y,'LineWidth',2)
                    xlabel('performance')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    mTitle = sprintf('performance %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(1));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(1));
                    matfigname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(1));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 1sd
                    figure
                    tp = getTPC(1,Tb_out_1sd);
                    histogram(tp,tpedge, 'Normalization', 'probability');
                    hold on
                    ylim([0 0.4])
                    %             pd = fitdist(tp', 'Normal');
                    %             x = xlim;
                    %             x = x(1):0.1:x(2);
                    %             y = pdf(pd,x);
                    %             plot(x,y,'LineWidth',2)
                    xlabel('performance')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    mTitle = sprintf('performance %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(2));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(2));
                    matfigname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(2));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 3 sd
                    figure
                    tp = getTPC(1,Tb_out_3sd);
                    histogram(tp, tpedge, 'Normalization', 'probability');
                    hold on
                    ylim([0 0.4])
                    %             pd = fitdist(tp', 'Normal');
                    %             x = xlim;
                    %             x = x(1):0.1:x(2);
                    %             y = pdf(pd,x);
                    %             plot(x,y,'LineWidth',2)
                    xlabel('performance')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    mTitle = sprintf('performance %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(3));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(3));
                    matfigname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(3));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 5sd
                    figure
                    tp = getTPC(1,Tb_out_5sd);
                    histogram(tp,tpedge, 'Normalization', 'probability');
                    hold on
                    ylim([0 0.4])
                    %             pd = fitdist(tp', 'Normal');
                    %             x = xlim;
                    %             x = x(1):0.1:x(2);
                    %             y = pdf(pd,x);
                    %             plot(x,y,'LineWidth',2)
                    xlabel('performance')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    mTitle = sprintf('performance %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(4));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(4));
                    matfigname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(4));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                    %% 10sd
                    figure
                    tp = getTPC(1,Tb_out_10sd);
                    histogram(tp, tpedge, 'Normalization', 'probability');
                    hold on
                    ylim([0 0.4])
                    %             pd = fitdist(tp', 'Normal');
                    %             x = xlim;
                    %             x = x(1):0.1:x(2);
                    %             y = pdf(pd,x);
                    %             plot(x,y,'LineWidth',2)
                    xlabel('performance')
                    ylabel('probability density function');
                    ax = gca; set(ax,'FontSize',12);
                    mTitle = sprintf('performance %1.3f; amp %1.3f; sd %1.3f', muPlus(j), seasonal_amp(k), sd(5));
                    tits = mtit(mTitle, 'fontsize', 12);
                    figname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus(j), seasonal_amp(k), sd(5));
                    matfigname = sprintf('performancePDF_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus(j), seasonal_amp(k), sd(5));
                    saveas(gcf, figname)
                    saveas(gcf, matfigname)
                end
            end
            
            %% plot Tb TPC frequency using spline interpolation of histogram
            if plotTbTPCspline
                plotTbTPC_spline(Tb_out_0sdNoise, muPlus(j), seasonal_amp(k), sd(1))
                plotTbTPC_spline(Tb_out_1sd, muPlus(j), seasonal_amp(k), sd(2))
                plotTbTPC_spline(Tb_out_3sd, muPlus(j), seasonal_amp(k), sd(3))
                plotTbTPC_spline(Tb_out_5sd, muPlus(j), seasonal_amp(k), sd(4))
                plotTbTPC_spline(Tb_out_10sd, muPlus(j), seasonal_amp(k), sd(5))
            end
        end
    end
end

end