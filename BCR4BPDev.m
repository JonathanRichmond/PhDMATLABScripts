%%% BCR4BPDev.jl
%%% Jonathan Richmond
%%% C: 19 February 2025
%%% U: 18 March 2025

clear
addpath('C:\Program Files\MATLAB\mice\src\mice\')
addpath('C:\Program Files\MATLAB\mice\lib\')
cspice_furnsh('naif0012.tls')
load("../PhDScripts/Output/BCR4BPDev.mat")
initialEpoch = cspice_et2utc(trajBCR4BPEEclipJ2000.t(1), 'C', 0);
MoonLoc = [MoonInitialState.x, MoonInitialState.y, MoonInitialState.z]./norm([MoonInitialState.x, MoonInitialState.y, MoonInitialState.z]);
SunLoc = [SunInitialState.x, SunInitialState.y, SunInitialState.z]./norm([SunInitialState.x, SunInitialState.y, SunInitialState.z]);

%% Body Data
% Earth
gmE = 3.9860043543609593E5; % Earth gravitational parameter [km^3/s^2]
mE = gmE/6.67384E-20; % Earth mass [kg]
RE = 6.3710083666666660E3; % Earth radius [km]

% Moon
gmm = 4.9028000661637961E3; % Moon gravitational parameter [km^3/s^2]
mm = gmm/6.67384E-20; % Moon mass [kg]
Rm = 1.7374000000000003E3; % Moon radius [km]

% Sun
gmS = 1.3271244004193930E11; % Sun gravitational parameter [km^3/s^2]
mS = gmS/6.67384E-20; % Sun mass [kg]
RS = 6.9600000000000000E5; % Sun radius [km]

% B1
gmB1 = gmE+gmm; %Barycenter gravitational parameter [km^3/s^2]
mB1 = gmB1/6.67384E-20; % Barycenter mass [kg]

%% Earth-Moon Data
mstarEM = (gmE+gmm)/6.67384E-20; % Characteristic mass [kg]
lstarEM = 3.8474799201129237E5; % Characteristic length [km]
tstarEM = sqrt(lstarEM^3/(gmE+gmm)); % Characteristic time [s]
muEM = gmm/(gmE+gmm); % Mass ratio

m_S = mS/mstarEM;
a_S = 1.4959789401764473E8/lstarEM;

MoonEclipJ2000 = rotToP1EclipJ2000(muEM, initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, trajCR3BPEclipJ2000.t, ones(length(trajCR3BPEclipJ2000.t), 1)*[1-muEM, 0, 0, 0, 0, 0]);

%% Sun-B1 Data
mstarSB1 = mstarEM*(m_S+1); % Characteristic mass [kg]
lstarSB1 = lstarEM*a_S; % Characteristic length [km]
tstarSB1 = sqrt(lstarSB1^3/(gmB1+gmS)); % Characteristic time [s]
muSB1 = mstarEM/mstarSB1;

a_M = (1-muEM)*lstarEM/lstarSB1;

%% Libration Points
g1 = muEM;
delg1 = 1;
while abs(delg1) > eps
    f = ((1-muEM)/((1-g1)^2))-(muEM/(g1^2))-1+muEM+g1;
    fprime = ((2*(1-muEM))/((1-g1)^3))+((2*muEM)/(g1^3))+1;
    g1New = g1-(f/fprime);
    delg1 = g1New-g1;
    g1 = g1New;
end
a1EM = 1-muEM-g1;

g2 = muEM;
delg2 = 1;
while abs(delg2) > eps
    f = ((1-muEM)/((1+g2)^2))+(muEM/(g2^2))-1+muEM-g2;
    fprime = ((-2*(1-muEM))/((1+g2)^3))-((2*muEM)/(g2^3))-1;
    g2New = g2-(f/fprime);
    delg2 = g2New-g2;
    g2 = g2New;
end
a2EM = 1-muEM+g2;

g3 = muEM;
delg3 = 1;
while abs(delg3) > eps
    f = ((1-muEM)/(g3^2))+(muEM/((1+g3)^2))-muEM-g3;
    fprime = ((-2*(1-muEM))/(g3^3))-((2*muEM)/((1+g3)^3))-1;
    g3New = g3-(f/fprime);
    delg3 = g3New-g3;
    g3 = g3New;
end
a3EM = -1*muEM-g3;

a45EM = 0.5-muEM;
b4 = sqrt(3)/2;
b5 = -b4;

%% Plots
fig1 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
set(Moon, 'DisplayName', "Moon")
scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "$L_{1}$")
scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "$L_{2}$")
scatter3(a3EM, 0, 0, 20, 'g', 'filled', 'd', 'DisplayName', "$L_{3}$")
scatter3(a45EM, b4, 0, 20, 'b', 'filled', 'd', 'DisplayName', "$L_{4}$")
scatter3(a45EM, b5, 0, 20, [1 0 1], 'filled', 'd', 'DisplayName', "$L_{5}$")
plot3(orbitCR3BP.x, orbitCR3BP.y, orbitCR3BP.z, 'DisplayName', "CR3BP Orbit")
scatter3(trajBCR4BPEM.x, trajBCR4BPEM.y, trajBCR4BPEM.z, 10*ones(length(trajBCR4BPEM.t), 1), angleColor(trajBCR4BPEM.theta4), 'filled', 'DisplayName', "BCR4BP EM Prop.")
plot3(validBCR4BPEM.x, validBCR4BPEM.y, validBCR4BPEM.z, 'DisplayName', "BCR4BP S$B_{1}$ Trans.")
plot3(validBCR4BPEM2.x, validBCR4BPEM2.y, validBCR4BPEM2.z, 'DisplayName', "ECLIPJ2000 Trans.")
axis equal
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
title("Earth-Moon: $\theta_{S,0}="+num2str(rad2deg(trajBCR4BPEM.theta4(1)), '%.0f')+"^{\circ}$", 'Interpreter', 'latex')
legend('Location', 'northeastoutside', 'Interpreter', 'latex')
phasemap;
phasebar('Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k')
view(3)
hold off
% exportgraphics(fig1, 'BCR4BPDev_1.png', 'BackgroundColor', 'k')

fig2 = figure("Position", [200 100 1200 750]);
hold on
B1 = plot3DBody("Earth", RE/lstarSB1, [1-muSB1, 0, 0]);
set(B1, 'DisplayName', "$B_{1}$")
fplot(@(t) a_M*sin(t)+1-muSB1, @(t) a_M*cos(t), 'w', 'DisplayName', "Lunar Orbit");
scatter3(trajBCR4BPSB1.x, trajBCR4BPSB1.y, trajBCR4BPSB1.z, 10*ones(length(trajBCR4BPSB1.t), 1), angleColor(trajBCR4BPSB1.theta2), 'filled', 'DisplayName', "BCR4BP EM Trans.")
plot3(validBCR4BPSB1.x, validBCR4BPSB1.y, validBCR4BPSB1.z, 'DisplayName', "BCR4BP S$B_{1}$ Prop.")
plot3(validBCR4BPSB12.x, validBCR4BPSB12.y, validBCR4BPSB12.z, 'DisplayName', "ECLIPJ2000 Trans.")
axis equal
grid on
xlabel("$\underline{x}$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
ylabel("$\underline{y}$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
zlabel("$\underline{z}$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
title("Sun-$B_{1}$: $\theta_{M,0}="+num2str(rad2deg(trajBCR4BPSB1.theta2(1)), '%.0f')+"^{\circ}$", 'Interpreter', 'latex')
legend('Location', 'northeastoutside', 'Interpreter', 'latex')
phasemap;
phasebar('Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k')
view(3)
hold off
% exportgraphics(fig2, 'BCR4BPDev_2.png', 'BackgroundColor', 'k')

fig3 = figure("Position", [200 100 1200 750]);
hold on
B1 = plot3DBody("Earth", RE/lstarEM*5, [0, 0, 0]);
set(B1, 'DisplayName', "$B_{1}$")
Moon = plot3DBody("Moon", Rm/lstarEM*10, MoonLoc);
set(Moon, 'DisplayName', "Moon")
Sun = plot3DBody("Sun", RS/lstarSB1*10, 2*SunLoc);
set(Sun, 'DisplayName', "Sun")
axis equal
grid on
xlabel("$X$ [meaningless]", 'Interpreter', 'latex')
ylabel("$Y$ [meaningless]", 'Interpreter', 'latex')
title("Ecliptic J2000 ($B_{1}$): "+initialEpoch, 'Interpreter', 'latex')
legend('Location', 'northeastoutside', 'Interpreter', 'latex')
set(gca, 'Color', 'k')
view(2)
hold off
% exportgraphics(fig3, 'BCR4BPDev_3.png', 'BackgroundColor', 'k')

fig4 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstarEM, [0, 0, 0]);
set(Earth, 'DisplayName', "Earth")
fplot(@(t) sin(t), @(t) cos(t), 'w', 'DisplayName', "Planar Lunar Orbit");
plot3(MoonEclipJ2000(:,1), MoonEclipJ2000(:,2), MoonEclipJ2000(:,3), 'w:', 'DisplayName', 'True Lunar Orbit')
scatter3(trajBCR4BPEEclipJ2000.x, trajBCR4BPEEclipJ2000.y, trajBCR4BPEEclipJ2000.z, 10*ones(length(trajBCR4BPEEclipJ2000.t), 1), angleColor(trajBCR4BPEEclipJ2000.theta4), 'filled', 'DisplayName', "BCR4BP EM Traj.")
plot3(trajBCR4BPSB1EEclipJ2000.x, trajBCR4BPSB1EEclipJ2000.y, trajBCR4BPSB1EEclipJ2000.z, 'DisplayName', "BCR4BP S$B_{1}$ Traj.")
plot3(trajCR3BPEclipJ2000.x, trajCR3BPEclipJ2000.y, trajCR3BPEclipJ2000.z, 'DisplayName', "CR3BP Traj.")
axis equal
grid on
xlabel("$X$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$Y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$Z$ [EM ndim]", 'Interpreter', 'latex')
title("EJ2000 (Earth): "+initialEpoch(1:11), 'Interpreter', 'latex')
legend('Location', 'northeastoutside', 'Interpreter', 'latex')
phasemap;
phasebar('Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k')
view(3)
% view(90, 0)
hold off
% exportgraphics(fig4, 'BCR4BPDev_4.png', 'BackgroundColor', 'k')

fig5 = figure("Position", [200 100 1200 750]);
hold on
Sun = plot3DBody("Sun", RS/lstarEM*5, [0, 0, 0]);
set(Sun, 'DisplayName', "Sun (x5)")
fplot(@(t) (lstarSB1/lstarEM)*sin(t), @(t) (lstarSB1/lstarEM)*cos(t), 'w', 'DisplayName', "$B_{1}$ Orbit");
scatter3(trajBCR4BPSEclipJ2000.x, trajBCR4BPSEclipJ2000.y, trajBCR4BPSEclipJ2000.z, 10*ones(length(trajBCR4BPSEclipJ2000.t), 1), angleColor(trajBCR4BPSEclipJ2000.theta4), 'filled', 'DisplayName', "BCR4BP EM Traj.")
axis equal
grid on
xlabel("$X$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$Y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$Z$ [EM ndim]", 'Interpreter', 'latex')
title("EJ2000 (Sun): "+initialEpoch(1:11), 'Interpreter', 'latex')
legend('Location', 'northeastoutside', 'Interpreter', 'latex')
phasemap;
phasebar('Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k')
view(3)
hold off
% exportgraphics(fig5, 'BCR4BPDev_5.png', 'BackgroundColor', 'k')

cspice_kclear