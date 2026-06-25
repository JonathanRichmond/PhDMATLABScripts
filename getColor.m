%%% getColor
%%% Jonathan LeFevre Richmond
%%% C: 24 June 2026

function color = getColor(colors, value, limits)

values = linspace(limits(1), limits(2), 1000);
i = 1;
while (i < 1000) && (values(i) < value)
    i = i+1;
end
i = min(i, 999);
v2 = values(i+1);
v1 = values(i);
c2 = colors(i+1,:);
c1 = colors(i,:);

color = max(0, min(1, c1+(((value-v1)*(c2-c1))/(v2-v1))));