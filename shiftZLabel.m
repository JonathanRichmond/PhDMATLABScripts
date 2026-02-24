%%% shiftZLabel
%%% Jonathan LeFevre Richmond
%%% C: 18 February 2026

function [] = shiftZLabel(ax, hz)

pos = get(hz, 'Position');
pos(1) = ax.XLim(1)-0.1*(ax.XLim(2)-ax.XLim(1));
pos(2) = ax.YLim(2)+0.1*(ax.YLim(2)-ax.YLim(1));
set(hz, 'Position', pos)