%%% Plot BCR4BP Trajectory
%%% Jonathan Richmond
%%% C: 11 June 2025
%%% U: 23 June 2025

clear
load('../PhDScripts/Output/EMSpatialOrbit.mat')

%% Earth-Moon Data
gmE = 3.9860043543609593E5; % Earth gravitational parameter [km^3/s^2]
mE = gmE/6.67384E-20; % Earth mass [kg]
RE = 6.3710083666666660E3; % Earth radius [km]

gmm = 4.9028000661637961E3; % Moon gravitational parameter [km^3/s^2]
mm = gmm/6.67384E-20; % Moon mass [kg]
Rm = 1.7374000000000003E3; % Moon radius [km]

mu = gmm/(gmE+gmm); % Mass ratio
mstar = (gmE+gmm)/6.67384E-20; % Characteristic mass [kg]
lstar = 3.8474799201129237E5; % Characteristic length [km]
tstar = sqrt(lstar^3/(gmE+gmm)); % Characteristic time [s]

%% Sun-B1 Data
gmS = 1.3271244004193930E11; % Sun gravitational parameter [km^3/s^2]
mS = gmS/6.67384E-20; % Sun mass [kg]
RS = 6.96E5; % Sun radius [km]

gmB1 = gmE+gmm; % B1 gravitational parameter [km^3/s^2]
mB1 = gmB1/6.67384E-20; % B1 mass [kg]

muSB1 = gmB1/(gmS+gmB1); % Mass ratio
mstarSB1 = (gmS+gmB1)/6.67384E-20; % Characteristic mass [kg]
lstarSB1 = 1.4959789401764473E8; % Characteristic length [km]
tstarSB1 = sqrt(lstarSB1^3/(gmS+gmB1)); % Characteristic time [s]

%% Plots
fig1 = figure("Position", [200 100 1200 750]);
hold on
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
plot3(BCR4BPOrbit.x, BCR4BPOrbit.y, BCR4BPOrbit.z, 'b', 'DisplayName', "BCR4BP 11:7 Syn.")
plot3(CR3BPCompOrbit.x, CR3BPCompOrbit.y, CR3BPCompOrbit.z, 'r', 'DisplayName', "CR3BP 11:7 Syn.")
axis equal
% axis
% axis([-0.16 1.88 -0.67 0.67])
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
title("Earth-Moon Rot.", 'Interpreter', 'latex')
legend('Location', 'bestoutside', 'Interpreter', 'latex');
% phasemap;
% phasebar('deg', 'Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k');
view(3)
hold off
% exportgraphics(fig1, 'EMSpatialOrbit_1.png', 'BackgroundColor', 'k')

fig2 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstarSB1, [1-muSB1, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moontheta = 0:(pi/180):2*pi;
Moonx = lstar.*cos(Moontheta);
Moony = lstar.*sin(Moontheta);
plot3(1-muSB1+Moonx./lstarSB1, Moony./lstarSB1, zeros(1, length(Moontheta)), 'w:', 'DisplayName', "Moon Orbit")
plot3(BCR4BPSB1Orbit.x, BCR4BPSB1Orbit.y, BCR4BPSB1Orbit.z, 'b', 'DisplayName', "BCR4BP 11:7 Syn.")
axis equal
% axis
% axis([-0.16 1.88 -0.67 0.67])
grid on
xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
zlabel("$z$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
title("Sun-$B_{1}$ Rot.", 'Interpreter', 'latex')
legend('Location', 'bestoutside', 'Interpreter', 'latex');
% phasemap;
% phasebar('deg', 'Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k');
view(3)
hold off
% exportgraphics(fig2, 'EMSpatialOrbit_2.png', 'BackgroundColor', 'k')