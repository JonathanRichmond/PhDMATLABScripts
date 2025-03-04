%%% plot2DProjections
%%% Jonathan Richmond
%%% C: 4 March 2025

function [] = plot2DProjections(ax, limits, vec1, vec2, vec3, darkness)

XL = ax.XLim;
YL = ax.YLim;
ZL = ax.ZLim;
plot3(vec1, vec2, ZL(limits(3))*ones(length(vec3)), 'Color', darkness*ones(1, 3), 'HandleVisibility', 'off')
plot3(XL(limits(1))*ones(length(vec1)), vec2, vec3, 'Color', darkness*ones(1, 3), 'HandleVisibility', 'off')
plot3(vec1, YL(limits(2))*ones(length(vec2)), vec3, 'Color', darkness*ones(1, 3), 'HandleVisibility', 'off')
axis([XL(1) XL(2) YL(1) YL(2) ZL(1) ZL(2)])