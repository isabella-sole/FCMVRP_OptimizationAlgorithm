function  [vertex_idxs, target_idxs, cs_idxs , h_figure] = create_map(CS, Target, xloc, yloc, target_lalbel, cs_lalbel)

vertex_idxs = [1:Target, Target+1:Target+CS]; % array of the vertexes' indexes 
target_idxs = [1:Target]; % array of the targets' indexes 
cs_idxs = [Target+1:Target+CS]; % array of the charging stations' indexes

% display the figure plotting the starting point, cs points and target points
%h = figure('visible','off');
h_figure = figure;
plot(xloc(1:Target), yloc(1:Target), 'rs', xloc(Target+1:length(xloc)), yloc(Target+1:length(yloc) ), 'bo');
text(xloc(1:Target),yloc(1:Target),target_lalbel,'VerticalAlignment','bottom','HorizontalAlignment','right');
text(xloc(Target+1:length(xloc)), yloc(Target+1:length(yloc)),cs_lalbel,'VerticalAlignment','bottom','HorizontalAlignment','right');
lgnd = legend('Target', 'Depot/Charging Station', 'Location', 'EastOutside');
lgnd.AutoUpdate = 'off';
xlim([0 15]);ylim([0 15])
title('Target and Charging Station Positions');
end