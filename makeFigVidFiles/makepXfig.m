function makepXfig(seasonal_amp, y, foodDensity, z, nnt, p_X, muPlus, nmap)

fid1 = figure;
figure(fid1)
hold on

for i = 1:length(muPlus)
plot(nnt, p_X(:,i), '-o', ...
    'color', nmap(i,1:3), ...
	'LineWidth', 2, ...
    'MarkerFaceColor' , nmap(i,1:3), ...
	'MarkerSize', 5)
end
title(strcat('+', int2str(muPlus), ' samp=', int2str(seasonal_amp(y)), ' X=', int2str(foodDensity(z))));
setpXProperties