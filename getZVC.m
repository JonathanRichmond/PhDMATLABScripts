%%% getZVC
%%% Jonathan Richmond
%%% C: 8 April 2025

function ZVC = getZVC(muEM, mS, aS, thetaS, H, axisLimits, numPoints, lineSpec)

[x, y] = meshgrid(linspace(axisLimits(1), axisLimits(2), numPoints), linspace(axisLimits(3), axisLimits(4), numPoints));
r_13 = sqrt((x+muEM).^2+y.^2);
r_23 = sqrt((x-1+muEM).^2+y.^2);
r_43 = sqrt((x-aS.*cos(thetaS)).^2+(y-aS.*sin(thetaS)).^2);
U = (1-muEM)./r_13+muEM./r_23+mS./r_43-mS.*(x.*cos(thetaS)+y.*sin(thetaS))./aS.^2+(1/2).*(x.^2+y.^2);
HEM = 2.*U;

[~, ZVC] = contour(x, y, HEM, [H, H], lineSpec);