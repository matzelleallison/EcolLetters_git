function makeTPCVid(seasonal_amp, y, foodDensity, z, nnt, T_b, muPlus, x, colors)
TPC = getTPC(1, -5:40);
frameRate = 60;

TP = colors(1:3);
pt = colors(4:6);
% [map,~,~] = brewermap(11,'*Spectral');

% % define colors for current and total TP
% muPlusof1_TP = map(1,:);    muPlusof1_pt = map(2,:);
% muPlusof2_TP = map(4,:);    muPlusof2_pt = map(5,:);
% muPlusof3_TP = map(8,:);    muPlusof3_pt = map(7,:);
% muPlusof4_TP = map(11,:);    muPlusof4_pt = map(10,:);


%% Plot final figure for testing before running video writer
% % ----- muPlus(1)
% fid1 = figure;
% figure(fid1)
% % plot all performance behind TPC curve
%     plot(T_b(:,1), getTPC(1,T_b(:,1)), 'o' , ...
%         'color' , muPlusof1_TP , ...
%         'MarkerFaceColor' , muPlusof1_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC,'k', 'LineWidth',5)
%     % plot one point on TPC curve
%     plot(T_b(1,1), getTPC(1,T_b(1,1)), 'o' , ...
%         'color', muPlusof1_pt, ...
%         'MarkerFaceColor' , muPlusof1_pt , ...
%         'MarkerSize', 15)
%     title(strcat('mu Tb = ', int2str(muPlus(1)), ' samp=', int2str(seasonal_amp(y))));
% %     xlabel('Temperature'); ylabel('Performance');
% %     ylim([0 1.1])
% %     ax = gca;
% %     box(ax,'off');
% %     % Set the remaining axes properties
% %     set(ax,'FontName','Gill Sans','FontSize',14);
%     setProperties
%     clear TPCvfname1 writerObj1 fid1 
% 
% % ----- muPlus(2)
% fid1 = figure;
% figure(fid1)
%     % plot all performance behind TPC curve
%     plot(T_b(:,2), getTPC(1,T_b(:,2)), 'o', ...
%         'color', muPlusof2_TP, ...
%         'MarkerFaceColor', muPlusof2_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC, 'k', 'LineWidth', 5)
%     % plot one point on TPC curve
%     plot(T_b(1,2), getTPC(1,T_b(1,2)), 'o' , ...
%         'color', muPlusof2_pt, ...
%         'MarkerFaceColor' , muPlusof2_pt , ...
%         'MarkerSize', 15)
%     title(strcat('mu Tb = ', int2str(muPlus(2)), ' samp=', int2str(seasonal_amp(y))));
%     setProperties
% clear TPCvfname1 writerObj1 fid1 
% 
% % ----- muPlus(3)
% fid1 = figure;
% figure(fid1)
% 
%     % plot all performance behind TPC curve
%     plot(T_b(:,3), getTPC(1,T_b(:,3)), 'o', ...
%         'color', muPlusof3_TP, ...
%         'MarkerFaceColor', muPlusof3_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC, 'k', 'LineWidth', 5)
%     % plot current point on TPC curve
%     plot(T_b(1,3), getTPC(1,T_b(1,3)), 'o' , ...
%         'color', muPlusof3_pt, ...
%         'MarkerFaceColor' , muPlusof3_pt , ...
%         'MarkerSize', 15)
%     title(strcat('mu Tb = ', int2str(muPlus(3)), ' samp=', int2str(seasonal_amp(y))));
%     setProperties
% clear TPCvfname1 writerObj1 fid1 
% 
% % ----- muPlus(4)
% fid1 = figure;
% figure(fid1)
% % plot previous performance behind TPC curve
%     plot(T_b(:,4), getTPC(1,T_b(:,4)), 'o', ...
%         'color', muPlusof4_TP, ...
%         'MarkerFaceColor', muPlusof4_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC, 'k', 'LineWidth', 5)
%     % plot current point on TPC curve
%     plot(T_b(1,4), getTPC(1,T_b(1,4)), 'o' , ...
%         'color', muPlusof4_pt, ...
%         'MarkerFaceColor' , muPlusof4_pt , ...
%         'MarkerSize', 15)
%     title(strcat('mu Tb = ', int2str(muPlus(4)), ' samp=', int2str(seasonal_amp(y))));
%     setProperties
% clear TPCvfname1 writerObj1 fid1 

%% Plot and save videos
% ----- muPlus(1)
TPCvfname1 = sprintf('TPC-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(x), seasonal_amp(y), foodDensity(z));
writerObj1 = VideoWriter(TPCvfname1); % Name it.
writerObj1.FrameRate = frameRate; % How many frames per second.
open(writerObj1);

fid1 = figure;
figure(fid1)

for b = 2:length(nnt)
    % plot previous performance behind TPC curve
    plot(T_b(b-1:b), getTPC(1,T_b(b-1:b)), 'o' , ...
        'color' , TP , ...
        'MarkerFaceColor' , TP, ...
        'MarkerSize', 20)
    hold on
    % plot TPC curve
    plot(-5:40, TPC,'k', 'LineWidth',5)
    % plot current pooint on TPC curve
    plot(T_b(b), getTPC(1,T_b(b)), 'o' , ...
        'color', pt, ...
        'MarkerFaceColor' , pt , ...
        'MarkerSize', 15)
    
    title(strcat('mu Tb = ', int2str(muPlus(x)), ' samp=', int2str(seasonal_amp(y))));
    setProperties
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj1, frame);
%     hold off
end
hold off
close(writerObj1); % Saves the movie.
clear TPCvfname1 writerObj1 fid1 

% % ----- muPlus(2)
% TPCvfname1 = sprintf('TPC-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(2), seasonal_amp(y), foodDensity(z));
% writerObj1 = VideoWriter(TPCvfname1); % Name it.
% writerObj1.FrameRate = frameRate; % How many frames per second.
% open(writerObj1);
% 
% fid1 = figure;
% figure(fid1)
% 
% for b = 1:length(nnt)
%     % plot previous performance behind TPC curve
%     plot(T_b(1:b,2), getTPC(1,T_b(1:b,2)), 'o', ...
%         'color', muPlusof2_TP, ...
%         'MarkerFaceColor', muPlusof2_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC, 'k', 'LineWidth', 5)
%     % plot current point on TPC curve
%     plot(T_b(b,2), getTPC(1,T_b(b,2)), 'o' , ...
%         'color', muPlusof2_pt, ...
%         'MarkerFaceColor' , muPlusof2_pt , ...
%         'MarkerSize', 15)
%     
%     title(strcat('mu Tb = ', int2str(muPlus(2)), ' samp=', int2str(seasonal_amp(y))));
%     setProperties
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj1, frame);
%     hold off
% end
% hold off
% close(writerObj1); % Saves the movie.
% clear TPCvfname1 writerObj1 fid1 
% 
% % ----- muPlus(3)
% TPCvfname1 = sprintf('TPC-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(3), seasonal_amp(y), foodDensity(z));
% writerObj1 = VideoWriter(TPCvfname1); % Name it.
% writerObj1.FrameRate = frameRate; % How many frames per second.
% open(writerObj1);
% 
% fid1 = figure;
% figure(fid1)
% 
% for b = 1:length(nnt)
%     % plot previous performance behind TPC curve
%     plot(T_b(1:b,3), getTPC(1,T_b(1:b,3)), 'o', ...
%         'color', muPlusof3_TP, ...
%         'MarkerFaceColor', muPlusof3_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC, 'k', 'LineWidth', 5)
%     % plot current point on TPC curve
%     plot(T_b(b,3), getTPC(1,T_b(b,3)), 'o' , ...
%         'color', muPlusof3_pt, ...
%         'MarkerFaceColor' , muPlusof3_pt , ...
%         'MarkerSize', 15)
%     
%     title(strcat('mu Tb = ', int2str(muPlus(3)), ' samp=', int2str(seasonal_amp(y))));
%     setProperties
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj1, frame);
%     hold off    
% end
% hold off
% close(writerObj1); % Saves the movie.
% clear TPCvfname1 writerObj1 fid1 
% 
% % ----- muPlus(4)
% TPCvfname1 = sprintf('TPC-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(4), seasonal_amp(y), foodDensity(z));
% writerObj1 = VideoWriter(TPCvfname1); % Name it.
% writerObj1.FrameRate = frameRate; % How many frames per second.
% open(writerObj1);
% 
% fid1 = figure;
% figure(fid1)
% 
% for b = 1:length(nnt)
% % plot previous performance behind TPC curve
%     plot(T_b(1:b,4), getTPC(1,T_b(1:b,4)), 'o', ...
%         'color', muPlusof4_TP, ...
%         'MarkerFaceColor', muPlusof4_TP, ...
%         'MarkerSize', 20)
%     hold on
%     % plot TPC curve
%     plot(-5:40, TPC, 'k', 'LineWidth', 5)
%     % plot current point on TPC curve
%     plot(T_b(b,4), getTPC(1,T_b(b,4)), 'o' , ...
%         'color', muPlusof4_pt, ...
%         'MarkerFaceColor' , muPlusof4_pt , ...
%         'MarkerSize', 15)
%     
%     title(strcat('mu Tb = ', int2str(muPlus(4)), ' samp=', int2str(seasonal_amp(y))));
%     setProperties
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj1, frame);
%     hold off
% end
% hold off
% close(writerObj1); % Saves the movie.
% clear TPCvfname1 writerObj1 fid1 