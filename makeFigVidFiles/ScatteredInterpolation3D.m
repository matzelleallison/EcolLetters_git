function ScatteredInterpolation3D(muPlus, seasonal_amp, foodDensity, sd, smpl_interval, Figs)
close all
doPlot = Figs(1); doSave = Figs(2);
sdGrp = {'noNoise', '1sdNoise', '3sdNoise', '5sdNoise', '10sdNoise'};

% s_a = []; T_b = []; s_d_Tb = []; mu_w = [];
for d = 1:length(foodDensity)
    fd = repmat(foodDensity(d), 4, 1);
    for a = 1:length(seasonal_amp)
        sa = repmat(seasonal_amp(a), 4, 1);
        s_a = []; T_b = []; s_d_Tb = []; mu_w = [];
        f_d = []; W_w = []; L_w = []; W_d = [];
        
        for b = 1:length(muPlus)
            Tb = repmat(muPlus(b), length(sdGrp), 1);
            
            for c = 1:length(sdGrp)
                
                %% T performance
                %             sd_Tb = repmat(sd(c), 4, 1);
                sd_Tb = sd(c);
                fname = sprintf('environment-env_plus_%1.3f-env_samp_%1.3f_%s.mat', ...
                    muPlus(b), seasonal_amp(a), sdGrp{c});
                load(fname)
                
                T = envSeason(:,1); % pull one env from n=10 replicate envs
                T_tp = getTPC(1, T);    % convert body temp to performance
%                 w(a).w(c, b) = mean(T_tp);     % mean long term performance
                muw(c) = mean(T_tp);
                s_d_Tb = [s_d_Tb; sd_Tb];
                
                %% DEB predictions
                
                fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f_%s.mat', ...
                    muPlus(b), seasonal_amp(a), foodDensity(d), sdGrp{c});
                load(fname)
                DEBdata = DEB_out(1);
                
                [Lw, Ww, Wd, pX, pM, pLw, pWw, Tb_out] = getRates(DEBdata, smpl_interval);
                
                maxLw(c) = max(Lw); maxWw(c) = max(Ww); maxWd(c) = max(Wd);
                
                
            end
            s_a = [s_a; sa];
            T_b = [T_b; Tb];
            mu_w = [mu_w; muw'];
            
            L_w = [L_w; maxLw'];
            W_w = [W_w; maxWw'];
            W_d = [W_d; maxWd'];
        end
        
        xlin = linspace(min(muPlus), max(muPlus), 50);
        wlin = linspace(min(sd), max(sd), 50);
        [X, Y] = meshgrid(xlin, wlin);
        
        %% thermal performance
        thermalPerformance = 1;
        if thermalPerformance
            f = scatteredInterpolant(T_b, s_d_Tb, mu_w);
            Z = f(X,Y);
            
            if doPlot
                % Tb,sd,TP 2D contour
                figure
                %     mesh(X,Y,Z) %interpolated
                surf(X,Y,Z,'EdgeColor','none') %interpolated
                axis tight;
                hold on
                colorbar
                caxis([0,1])
                view(0,90)
                colormap Jet
                
                xlabel('mean temperature');
                ylabel('standard deviation in temperature');
                title(strcat('seasonal amplitude = ', int2str(seasonal_amp(a))));
                ax = gca; set(ax,'FontSize',12);

                if doSave
                    figname = sprintf('mean_Tb-std_Tb-w_samp-%1.3f.tiff', seasonal_amp(a));     
                    matfigname = sprintf('mean_Tb-std_Tb-w_samp-%1.3f.fig', seasonal_amp(a));   
                    saveas(gcf, figname)
                    saveas(gcf,matfigname)
                end
            end
        end
        
        %% DEB outs
       
        
        % shell length
        sI_Lw = scatteredInterpolant(T_b, s_d_Tb, L_w);
        Z_Lw = sI_Lw(X, Y);
        if doPlot
            figure
            surf(X, Y, log(Z_Lw),'EdgeColor','none')
            axis tight;
            hold on
            h = colorbar;
            title(h, 'Shell length','FontSize',10)
            caxis([1.8,3.3])   % max length (f = 0.99, T = Topt) = 32.40; (f =  0.5 , T = Topt) = 16.89; (f =  0.3 , T = Topt) = 10.22
            view(0,90)
            colormap Jet
            xlabel('mean temperature','FontSize',10);
            ylabel('standard deviation in temperature','FontSize',10);
            title(strcat('seasonal amplitude = ', int2str(seasonal_amp(a)), '; food density = ', int2str(foodDensity(d))),'FontSize',10);
            ax = gca; set(ax,'FontSize',12);
            
%             fig = gcf;
%             fig.Units = 'inches';
%             fig.Position = [0 0 2.5 2.5];
%             fig.PaperUnits = 'inches';
%             fig.PaperPosition = [0 0 2.5 2.5];
%             
%             ax = gca;
%             ax.FontSize = 12;
            %             outerpos = ax.OuterPosition;
            %             ti = ax.TightInset;
            %             left = outerpos(1) + ti(1);
            %             bottom = outerpos(2) + ti(2);
            %             ax_width = outerpos(3) - ti(1) - ti(3);
            %             ax_height = outerpos(4) - ti(2) - ti(4);
            %             ax.Position = [left bottom ax_width ax_height];
            if doSave
                figname = sprintf('mean_Tb-std_Tb-Lw_samp-%1.3f_food-%1.3f.tiff', seasonal_amp(a), foodDensity(d)); 
                matfigname = sprintf('mean_Tb-std_Tb-Lw_samp-%1.3f_food-%1.3f.fig', seasonal_amp(a), foodDensity(d));   
                
                saveas(gcf, figname)
                saveas(gcf,matfigname)
            end
        end
            
        % wet weight
        sI_Ww = scatteredInterpolant(T_b, s_d_Tb, W_w);
            Z_Ww = sI_Ww(X, Y);
            if doPlot
                figure
                surf(X, Y, log(Z_Ww),'EdgeColor','none')
                axis tight;
                hold on
                colorbar
                h = colorbar;
                title(h, 'Wet weight')
                caxis([1.4,6.2])  % max wet weight (f = 0.99, T = Topt) = 704.20; (f =  0.5 , T = Topt) = 96.64; (f =  0.3 , T = Topt) = 21.17
                view(0,90)
                colormap Jet
                xlabel('mean temperature');
                ylabel('standard deviation in temperature');
                title(strcat('seasonal amplitude = ', int2str(seasonal_amp(a)), '; food density = ', int2str(foodDensity(d))));
%                 set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
                ax = gca; set(ax,'FontSize',12);
                
                if doSave
                    figname = sprintf('mean_Tb-std_Tb-Ww_samp-%1.3f_food-%1.3f.tiff', seasonal_amp(a), foodDensity(d));
                    matfigname = sprintf('mean_Tb-std_Tb-Ww_samp-%1.3f_food-%1.3f.fig', seasonal_amp(a), foodDensity(d));
                    
                    saveas(gcf, figname)
                    saveas(gcf,matfigname)
                end
            end
        
            % dry weight
            sI_Wd = scatteredInterpolant(T_b, s_d_Tb, W_d);
            Z_Wd = sI_Wd(X, Y);
            if doPlot
                figure
                surf(X, Y, log(Z_Wd),'EdgeColor','none')
                axis tight;
                hold on
                colorbar
                h = colorbar;
                title(h, 'Dry weight')
                caxis([0,4.2])   % max dry weight (f = 0.99, T = Topt) = 105.11; (f =  0.5 , T = Topt) = 11.74; (f = 0.3, T=Topt)=2.32
                view(0,90)
                colormap Jet
                xlabel('mean temperature');
                ylabel('standard deviation in temperature');
                title(strcat('seasonal amplitude = ', int2str(seasonal_amp(a)), '; food density = ', int2str(foodDensity(d))));
                ax = gca; set(ax,'FontSize',12);
                
                if doSave
                    figname = sprintf('mean_Tb-std_Tb-Wd_samp-%1.3f_food-%1.3f.tiff', seasonal_amp(a), foodDensity(d));
                    matfigname = sprintf('mean_Tb-std_Tb-Wd_samp-%1.3f_food-%1.3f.fig', seasonal_amp(a), foodDensity(d));
                    
                    saveas(gcf, figname)
                    saveas(gcf,matfigname)
                end
            end
        
    end
end
end




















%% wet weight
% % low food density
% for y = 1:length(seasonal_amp)
%     for x = 1:length(muPlus)
%         fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', muPlus(x), seasonal_amp(y), foodDensity(1));
%         load(fname)
%
%         % pull out Ww data from each replicate model run
%         for i = 1:10
%             data = DEB_out(i); mainTitle = '';
%             Ww_wiTemp(i) = real(max(data.W_w)); % store in matrix
%         end
%         meanWw_wiTemp(x,y) = mean(Ww_wiTemp);    % store in m,n matrix where inc avg temp rowwise, inc amplitude column-wise
%     end
% end
% surf(X, Y, meanWw_wiTemp)
%
%
% V = meanWw_wiTemp;
% V = reshape(V,[],1);
%
% [xq,yq] = meshgrid(linspace(0,40,200), linspace(0,15,200));
% vq = griddata(X,Y,V,xq,yq,'cubic');
% figure
% mesh(xq,yq,vq);
% hold on
% plot3(X,Y,V,'o','MarkerFaceColor' , 'auto', 'MarkerSize', 15);
%
% xlabel('Average body temperature')
% ylabel('seasonal amplitude')
% zlabel('max wet weight')
% title('Low food')
%
% % med food density
% for y = 1:length(seasonal_amp)
%     for x = 1:length(muPlus)
%         fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', muPlus(x), seasonal_amp(y), foodDensity(2));
%         load(fname)
%
%         % pull out Ww data from each replicate model run
%         for i = 1:10
%             data = DEB_out(i); mainTitle = '';
%             Ww_wiTemp(i) = real(max(data.W_w)); % store in matrix
%         end
%         meanWw_wiTemp(x,y) = mean(Ww_wiTemp);    % store in m,n matrix where inc avg temp rowwise, inc amplitude column-wise
%     end
% end
% V = meanWw_wiTemp;
% V = reshape(V,[],1);
%
% [xq,yq] = meshgrid(linspace(0,40,200), linspace(0,15,200));
% vq = griddata(X,Y,V,xq,yq,'cubic');
% figure
% mesh(xq,yq,vq);
% hold on
% plot3(X,Y,V,'o','MarkerFaceColor' , 'auto', 'MarkerSize', 15);
%
% xlabel('Average body temperature')
% ylabel('seasonal amplitude')
% zlabel('max wet weight')
% title('Med food')
%
% % high food density
% for y = 1:length(seasonal_amp)
%     for x = 1:length(muPlus)
%         fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', muPlus(x), seasonal_amp(y), foodDensity(3));
%         load(fname)
%
%         % pull out Ww data from each replicate model run
%         for i = 1:10
%             data = DEB_out(i); mainTitle = '';
%             Ww_wiTemp(i) = real(max(data.W_w)); % store in matrix
%         end
%         meanWw_wiTemp(x,y) = mean(Ww_wiTemp);    % store in m,n matrix where inc avg temp rowwise, inc amplitude column-wise
%     end
% end
% V = meanWw_wiTemp;
% V = reshape(V,[],1);
%
% [xq,yq] = meshgrid(linspace(0,40,200), linspace(0,15,200));
% vq = griddata(X,Y,V,xq,yq,'cubic');
% figure
% mesh(xq,yq,vq);
% hold on
% plot3(X,Y,V,'o','MarkerFaceColor' , 'auto', 'MarkerSize', 15);
%
% xlabel('Average body temperature')
% ylabel('seasonal amplitude')
% zlabel('max wet weight')
% title('High food')

