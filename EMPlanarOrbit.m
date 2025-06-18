%%% Plot BCR4BP Trajectory
%%% Jonathan Richmond
%%% C: 4 June 2025
%%% U: 9 June 2025

clear
load('../PhDScripts/Output/EMPlanarOrbit.mat')

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

%% Plot
fig1 = figure("Position", [200 100 1200 750]);
hold on
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
plot3(BCR4BPOrbit1.x, BCR4BPOrbit1.y, BCR4BPOrbit1.z, 'b', 'DisplayName', "BCR4BP 1:1(a) Syn.")
plot3(CR3BPGuessOrbit1.x, CR3BPGuessOrbit1.y, CR3BPGuessOrbit1.z, 'g', 'DisplayName', "CR3BP 1:1 Sid.")
plot3(CR3BPCompOrbit1.x, CR3BPCompOrbit1.y, CR3BPCompOrbit1.z, 'r', 'DisplayName', "CR3BP 1:1 Syn.")
plot3(BCR4BPOrbit2.x, BCR4BPOrbit2.y, BCR4BPOrbit2.z, 'c', 'DisplayName', "BCR4BP 1:1(b) Syn.")
plot3(CR3BPGuessOrbit2.x, CR3BPGuessOrbit2.y, CR3BPGuessOrbit2.z, 'g', 'DisplayName', "CR3BP 1:1 Sid.")
plot3(CR3BPCompOrbit2.x, CR3BPCompOrbit2.y, CR3BPCompOrbit2.z, 'r', 'DisplayName', "CR3BP 1:1 Syn.")
axis equal
% axis([-0.16 1.88 -0.67 0.67])
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
title("Earth-Moon $L_{1}$ Lyapunovs", 'Interpreter', 'latex')
legend('Location', 'bestoutside', 'Interpreter', 'latex');
% phasemap;
% phasebar('deg', 'Location', 'northeast', 'Size', 0.275)
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'EMPlanarOrbit_1.png', 'BackgroundColor', 'k')