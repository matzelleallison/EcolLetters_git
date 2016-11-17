function plot_BGYRhist(xi, f, i, Xlabel, Ylabel, Title)

colorOrder = [0 0.447 0.741;
    0.466 0.674 0.188;
    0.929 0.694 0.125;
    0.850 0.325 0.098];

legendval = {'mu = 5', 'mu = 15', 'mu = 25', 'mu = 35'};

h1 = area(xi,f,'FaceAlpha', 0.3, 'FaceColor', colorOrder(i,:));
hold on
h2 = plot(xi,f,'LineWidth', 2, 'color', colorOrder(i,:));

% legend(h2,legendval(i));

xlabel(Xlabel)
ylabel(Ylabel)
title(Title)


