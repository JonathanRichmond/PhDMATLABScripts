%%% PeriapsisMap.jl
%%% Jonathan Richmond
%%% C: 24 October 2025

clear
load('../PhDScripts/Output/PeriapsisMap.mat')

%% Earth-Moon Data
gmE = 3.9860043543609593E5; % Earth gravitational parameter [km^3/s^2]
mE = gmE/6.67384E-20; % Earth mass [kg]
RE = 6.371008366666666E3; % Earth radius [km]

gmm = 4.9028000661637961E3; % Moon gravitational parameter [km^3/s^2]
mm = gmm/6.67384E-20; % Moon mass [kg]
Rm = 1.7374000000000003E3; % Moon radius [km]

mu = gmm/(gmE+gmm); % Mass ratio
mstar = (gmE+gmm)/6.67384E-20; % Characteristic mass [kg]
lstar = 3.8474799201129237E5; % Characteristic length [km]
tstar = sqrt(lstar^3/(gmE+gmm)); % Characteristic time [s]

RH = 1.496557260502363E6; % Hills radius [km]

%% Sun-B1 Data
gmS = 1.327124400419393E11; % Sun gravitational parameter [km^3/s^2]
mS = gmS/6.67384E-20; % Sun mass [kg]
RS = 6.96E5; % Sun radius [km]

gmB1 = gmE+gmm; % B1 gravitational parameter [km^3/s^2]
mB1 = gmB1/6.67384E-20; % B1 mass [kg]

muSB1 = gmB1/(gmS+gmB1); % Mass ratio
mstarSB1 = (gmS+gmB1)/6.67384E-20; % Characteristic mass [kg]
lstarSB1 = 1.4959789401764473E8; % Characteristic length [km]
tstarSB1 = sqrt(lstarSB1^3/(gmS+gmB1)); % Characteristic time [s]

%% Sun-B1 Periapsis Map
fig1 = figure("Position", [200 100 1200 750]);
hold on
colors = zeros(length(flags), 3);
x = zeros(length(flags), 1);
y = zeros(length(flags), 1);
parfor j = 1:length(points)
    switch flags(j)
        case 0
            colors(j,:) = [0, 0, 0];
        case {1, 2}
            colors(j,:) = [1, 0, 0];
        case 3
            colors(j,:) = 0.7*[1, 1, 1];
        case 4
            colors(j,:) = [0, 0, 1];
        case 5
            colors(j,:) = [0, 1, 1];
        case 6
            colors(j,:) = [0, 1, 0];
    end
    x(j) = points{j}(1);
    y(j) = points{j}(2);
end
scatter(x, y, 3, colors, 'filled', 'HandleVisibility', 'off')
B1 = plot3DBody("Earth", RE/lstarSB1, [1-muSB1, 0, 0]);
set(B1, 'DisplayName', "Earth")
scatter(lstar/lstarSB1*cos(moonAngle)+1-muSB1, lstar/lstarSB1*sin(moonAngle), 20, 'w', 'filled', 'DisplayName', "Moon")
plot(lstar/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, lstar/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w:', 'DisplayName', "Moon Radius")
% plot(RH/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
scatter(nan, nan, 20, [0, 0, 1], 'filled', 'DisplayName', "Direct Escape")
scatter(nan, nan, 20, [0, 1, 1], 'filled', 'DisplayName', "Escape + 1 Rev")
scatter(nan, nan, 20, [0, 1, 0], 'filled', 'DisplayName', "Escape + 2 Revs")
axis equal
grid on
xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
title("Periapsis Map $JC_{EM}="+num2str(JC)+"$ $|$ $\theta_{M}="+num2str(moonAngle*180/pi)+"^\circ$", 'Interpreter', 'latex')
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'PeriapsisMap_1.png', 'BackgroundColor', 'k')

%% Trajectory Plots
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% B1 = plot3DBody("Earth", RE/lstarSB1, [1-muSB1, 0, 0]);
% set(B1, 'DisplayName', "Earth")
% plot(lstar/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, lstar/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w:', 'DisplayName', "Moon Radius")
% plot(RH/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% plot(Traj41.x, Traj41.y, 'g', 'DisplayName', "Trajectory")
% scatter(Traj41.x(1), Traj41.y(1), 50, 'w', '*', 'DisplayName', "S$B_{1}$ Peri.")
% axis equal
% grid on
% xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
% title("Sun-$B_{1}$ Rot.", 'Interpreter', 'latex')
% leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig2, 'PeriapsisMap_2.png', 'BackgroundColor', 'k')
% 
% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% plot(Traj12.x, Traj12.y, 'g', 'DisplayName', "Trajectory")
% scatter(Traj12.x(1), Traj12.y(1), 50, 'w', '*', 'DisplayName', "S$B_{1}$ Peri.")
% axis equal
% grid on
% xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig3, 'PeriapsisMap_3.png', 'BackgroundColor', 'k')