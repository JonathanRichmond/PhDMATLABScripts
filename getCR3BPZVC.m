%%% getCR3BPZVC
%%% Jonathan Richmond
%%% C: 17 June 2025

function ZVC = getCR3BPZVC(mu, JC, axisLimits, numPoints, lineSpec)

[x, y] = meshgrid(linspace(axisLimits(1), axisLimits(2), numPoints), linspace(axisLimits(3), axisLimits(4), numPoints));
r_13 = sqrt((x+mu).^2+y.^2);
r_23 = sqrt((x-1+mu).^2+y.^2);
U = (1-mu)./r_13+mu./r_23+(1/2).*(x.^2+y.^2);
H = 2.*U;

[~, ZVC] = contour(x, y, H, [JC, JC], lineSpec);