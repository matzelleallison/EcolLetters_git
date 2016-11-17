set(0,'defaultfigurecolor',[1 1 1])    
xlabel('Temperature'); ylabel('Performance');
    ylim([0 1.1])
    ax = gca;
    box(ax,'off');
    % Set the remaining axes properties
    set(ax,'FontName','Gill Sans','FontSize',14);