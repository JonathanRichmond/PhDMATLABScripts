%%% shiftZLabel
%%% Jonathan LeFevre Richmond
%%% C: 18 February 2026
%%% U: 13 May 2026

function [] = shiftZLabel(ax, hz, shift)

pos = get(hz, 'Position');
pos(1) = ax.XLim(1)-shift*(ax.XLim(2)-ax.XLim(1));
pos(2) = ax.YLim(2)+shift*(ax.YLim(2)-ax.YLim(1));
set(hz, 'Position', pos)