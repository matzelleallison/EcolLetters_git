function makeffig(X)
col = [0.254901960784314,0.670588235294118,0.364705882352941];
    
fid1 = figure;
figure(fid1)
hold on

plot(0:150, get_f(0:150, 6.12), 'k', 'LineWidth', 5)

for i = 1:length(X)
    plot(X(i), get_f(X(i), 6.12), 'o', ...
        'color', col, ...
        'MarkerFaceColor', col, ...
        'MarkerSize', 15)
end
title(int2str(X))

xlabel('Food Density'); ylabel('Functional response');
ax = gca;
box(ax, 'off');
set(ax, 'FontName', 'Gill Sans', 'FontSize', 14);