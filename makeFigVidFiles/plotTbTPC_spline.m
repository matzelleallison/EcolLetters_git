function plotTbTPC_spline(Tbdata, muPlus, seasonal_amp, sd)

figure
%     histogram(Tbdata, 'EdgeColor', 'none')
%     hist(Tbdata)
%     h = findobj(gca,'Type','patch');
% 
%     h.FaceColor = [0 0.6 1];
%     h.EdgeColor = 'w';
    hold on
    [heights, centers] = hist(Tbdata);
        c1 = centers(1)-(centers(2)-centers(1));    c2 = centers(end)+(centers(2)-centers(1));
    ncenters = [c1 centers c2];
    nheights = [0 heights 0];
    spline1 = csape(ncenters,[0,nheights,0],'complete'); 
%     fnplt(spline1,'-',2)
    points = fnplt(spline1,'-',2);
    
    z = fnzeros(spline1);
    
    X = points(1,:);
    Y = points(2,:);

    if size(z) == 2
        plot(X,Y, 'LineWidth', 2)
    else
        a = max(find(z(2,:) < muPlus));
        b = min(find(z(2,:) > muPlus));
        if isempty(b)
            b = max(points(1,:));
        else
            b = z(2,b);
        end
        a = z(2,a);
        
        idx = find(X>a & X<b);
        idx = [idx(1)-1 idx idx(end)+1];

%         xymin = min(find(points(1,:) == a));
%         xymax = max(find(points(1,:) == b));

        plot(X(idx), Y(idx), 'LineWidth', 2)
    end
    
    xlim([-10 45])
    ylim([0 length(Tbdata)])
    
    xlabel('body temperature')
    ylabel('frequency');
    ax = gca; set(ax,'FontSize',12);

    mTitle = sprintf('Tb %1.3f; amp %1.3f; sd %1.3f', muPlus, seasonal_amp, sd);
    tits = mtit(mTitle, 'fontsize', 12);
    figname = sprintf('TbFreq_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus, seasonal_amp, sd);
    matfigname = sprintf('TbFreq_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus, seasonal_amp, sd);
    
        saveas(gcf, figname)
        saveas(gcf, matfigname)

figure
    tp = getTPC(1,Tbdata);
%     histogram(tp, 'NumBins', 100, 'EdgeColor', 'none')
% hist(tp*100)
    hold on
    [heights, centers] = hist(tp*100);
%         c1 = centers(1)-(centers(2)-centers(1));    c2 = centers(end)+(centers(2)-centers(1));
%     ncenters = [centers(1)-(centers(2)-centers(1)) centers centers(end)+(centers(2)-centers(1))];
%     nheights = [0 heights 0];
    spline1 = csape(centers,[0,heights,0],'complete'); 
    
    
     points = fnplt(spline1,'-',2);
    
    z = fnzeros(spline1);
    
    X = points(1,:);
    Y = points(2,:);

    if size(z) == 2
        plot(X,Y, 'LineWidth', 2)
    else
        a = max(find(z(2,:) < nanmean(tp*100)));
        b = min(find(z(2,:) > nanmean(tp*100)));
        if isempty(b)
            b = max(points(1,:));
        else
            b = z(2,b);
        end
        if isempty(a)
            a = min(points(1,:));
        else
            a = z(2,a);
        end
        
        idx = find(X>a & X<b);
        idx = [idx(1)-1 idx idx(end)+1];

%         xymin = min(find(points(1,:) == a));
%         xymax = max(find(points(1,:) == b));

        plot(X(idx), Y(idx), 'LineWidth', 2)
    end
%     fnplt(spline1,'-',2)
    xlim([0 100])
    ylim([0 length(Tbdata)])
    
    xlabel('performance')
    ylabel('frequency');
    ax = gca; set(ax,'FontSize',12);
            
    mTitle = sprintf('performance Tb %1.3f; amp %1.3f; sd %1.3f', muPlus, seasonal_amp, sd);
    tits = mtit(mTitle, 'fontsize', 12);
    figname = sprintf('performanceFreq_Tb_%1.3f-amp_%1.3f-sd_%1.3f.tiff', muPlus, seasonal_amp, sd);
    matfigname = sprintf('performanceFreq_Tb_%1.3f-amp_%1.3f-sd_%1.3f.fig', muPlus, seasonal_amp, sd);
    
        saveas(gcf, figname)
        saveas(gcf, matfigname)