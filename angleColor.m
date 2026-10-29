%%% angleColor
%%% Jonathan LeFevre Richmond
%%% C: 19 February 2025
%%% U: 28 October 2026

function [colors] = angleColor(values)

% valuesWrapped = phasewrap(values);
colorMap = phasemap(1001);
valueMap = linspace(0, 2*pi, 1001);
colors = zeros(length(values), 3);
for j = 1:length(values)
    % [~, colorIndex] = min(abs(valueMap-valuesWrapped(j)));
    [~, colorIndex] = min(abs(valueMap-values(j)));
    colors(j,:) = colorMap(colorIndex,:);
end