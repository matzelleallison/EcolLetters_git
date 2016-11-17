function makefoodxWwfig(seasonal_amp, y, foodDensity, z, nnt, W_w, muPlus, nmap)

fid1 = figure;
figure(fid1)
hold on

for i = 1:length(muPlus)
plot(nnt, W_w(:,i), ':', ...
    'color', nmap(i,1:3), ...
	'LineWidth', 2, ...
    'MarkerFaceColor' , nmap(i,1:3), ...
	'MarkerSize', 2)
end
ylim([0 30])
title(strcat('+', int2str(muPlus), ' samp=', int2str(seasonal_amp(y)), ' X=', int2str(foodDensity(z))));
setWwProperties