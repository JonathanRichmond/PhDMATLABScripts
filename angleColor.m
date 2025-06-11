%%% angleColor
%%% Jonathan Richmond
%%% C: 19 February 2025

function [colors] = angleColor(values)

valuesWrapped = phasewrap(values);
colorMap = phasemap(1001);
valueMap = linspace(-pi, pi, 1001);
colors = zeros(length(values), 3);
for j = 1:length(values)
    [~, colorIndex] = min(abs(valueMap-valuesWrapped(j)));
    colors(j,:) = colorMap(colorIndex,:);
end