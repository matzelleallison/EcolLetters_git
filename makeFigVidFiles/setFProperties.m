set(0,'defaultfigurecolor',[1 1 1])
% setTbProperties
xlabel('Days'); ylabel('Fecundity');
%     ylim([-5 50])
%     xlim([0 length(nnt)])
    ax = gca;
    box(ax,'off');
    % Set the remaining axes properties
    set(ax,'FontName','Gill Sans','FontSize',14);