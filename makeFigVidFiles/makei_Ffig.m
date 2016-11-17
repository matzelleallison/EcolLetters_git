function makei_Ffig(seasonal_amp, y, foodDensity, z, nt, i_F, muPlus, nmap)

i_F(i_F == 0) = NaN;

fid1 = figure;
figure(fid1)
hold on

tF = linspace(1,nt-365,10);

for i = 1:length(muPlus)
plot(tF, i_F(:,i), 'o', ...
    'color', nmap(i,1:3), ...
	'LineWidth', 2, ...
    'MarkerFaceColor' , nmap(i,1:3), ...
	'MarkerSize', 8)
end
ylim([0 4.5E7])
title(strcat('+', int2str(muPlus), ' samp=', int2str(seasonal_amp(y)), ' X=', int2str(foodDensity(z))));
setFProperties