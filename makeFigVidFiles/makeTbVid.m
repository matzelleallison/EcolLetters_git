function makeTbVid(seasonal_amp, y, foodDensity, z, nnt, T_b, muPlus, x, colors)
% colors: 2 by 4 matrix with 'TP' and 'pt' color
frameRate = 60;

TP = colors(1:3);
pt = colors(4:6);

%% Plot and save videos
% ----- muPlus(1)
Tbvfname = sprintf('Tb-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(x), seasonal_amp(y), foodDensity(z));
writerObj = VideoWriter(Tbvfname); % Name it.
writerObj.FrameRate = frameRate; % How many frames per second.
open(writerObj);

fid = figure;
figure(fid)
hold on

for b = 2:length(nnt)
    plot(nnt(b-1:b), T_b(b-1:b), '-' , ...
        'color' , TP , ...
        'MarkerFaceColor' , TP, ...
        'LineWidth', 2)
    plot(nnt(b), T_b(b), 'o' , ...
        'color' , pt , ...
        'MarkerFaceColor' , pt, ...
        'MarkerSize', 7)
    hold on
    title(strcat('Tb = ', int2str(muPlus(x)), ' samp=', int2str(seasonal_amp(y))));
    setTbProperties
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end
hold off
close(writerObj); % Saves the movie.
clear TPCvfname writerObj fid1

% % ----- muPlus(2)
% Tbvfname = sprintf('Tb-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(2), seasonal_amp(y), foodDensity(z));
% writerObj = VideoWriter(Tbvfname); % Name it.
% writerObj.FrameRate = frameRate; % How many frames per second.
% open(writerObj);
% 
% fid = figure;
% figure(fid)
% hold on
% for b = 1:length(nnt)
%     plot(nnt(1:b), T_b(1:b,2), 'o-' , ...
%         'color' , muPlusof2_TP , ...
%         'MarkerFaceColor' , muPlusof2_TP, ...
%         'MarkerSize', 5, ...
%         'LineWidth', 2)
%     hold on
%     title(strcat('Tb = ', int2str(muPlus(2)), ' samp=', int2str(seasonal_amp(y))));
%     setTbProperties
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj, frame);
% end
% hold off
% close(writerObj); % Saves the movie.
% clear TPCvfname writerObj fid1
% 
% % ----- muPlus(3)
% Tbvfname = sprintf('Tb-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(3), seasonal_amp(y), foodDensity(z));
% writerObj = VideoWriter(Tbvfname); % Name it.
% writerObj.FrameRate = frameRate; % How many frames per second.
% open(writerObj);
% 
% fid = figure;
% figure(fid)
% hold on
% for b = 1:length(nnt)
%     plot(nnt(1:b), T_b(1:b,3), 'o-' , ...
%         'color' , muPlusof3_TP , ...
%         'MarkerFaceColor' , muPlusof3_TP, ...
%         'MarkerSize', 5, ...
%         'LineWidth', 2)
%     hold on
%     title(strcat('Tb = ', int2str(muPlus(3)), ' samp=', int2str(seasonal_amp(y))));
%     setTbProperties
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj, frame);
% end
% hold off
% close(writerObj); % Saves the movie.
% clear TPCvfname writerObj fid1
% 
% % ----- muPlus(4)
% Tbvfname = sprintf('Tb-envplus_%1.3f-env_samp_%1.3f-env_X_%1.3f.avi', muPlus(4), seasonal_amp(y), foodDensity(z));
% writerObj = VideoWriter(Tbvfname); % Name it.
% writerObj.FrameRate = frameRate; % How many frames per second.
% open(writerObj);
% 
% fid = figure;
% figure(fid)
% hold on
% for b = 1:length(nnt)
%     plot(nnt(1:b), T_b(1:b,4), 'o-' , ...
%         'color' , muPlusof4_TP , ...
%         'MarkerFaceColor' , muPlusof4_TP, ...
%         'MarkerSize', 5, ...
%         'LineWidth', 2)
%     hold on
%     title(strcat('Tb = ', int2str(muPlus(4)), ' samp=', int2str(seasonal_amp(y))));
%     setTbProperties
%     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%     writeVideo(writerObj, frame);
% end
% hold off
% close(writerObj); % Saves the movie.
% clear TPCvfname writerObj fid1