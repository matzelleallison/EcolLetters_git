function makeTbfig(seasonal_amp, y, foodDensity, z, nnt, T_b, muPlus, nmap)

fid1 = figure;
figure(fid1)
hold on

for i = 1:length(muPlus)
plot(nnt, T_b(:,i), '-', ...
    'color', nmap(i,1:3), ...
	'LineWidth', 2) ...
%     'MarkerFaceColor' , nmap(i,1:3), ...
% 	'MarkerSize', 5)
end
title(strcat('+', int2str(muPlus), ' samp=', int2str(seasonal_amp(y)), ' X=', int2str(foodDensity(z))));
setTbProperties