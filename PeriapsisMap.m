%%% PeriapsisMap.jl
%%% Jonathan Richmond
%%% C: 24 October 2025
%%% U: 8 November 2025

clear
load('../PhDScripts/Output/PeriapsisMap.mat')
n = length(flagsEM);

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

%% Color Map
colorMap = viridis(7);
colorMap(8,:) = [0, 0, 0];
colorMap(9,:) = [1, 0, 0];
colorMap(10,:) = 0.7*[1, 1, 1];

%% Earth-Moon Periapsis Map
fig1 = figure("Position", [200 100 1200 750]);
hold on
x = zeros(n, 1);
y = zeros(n, 1);
parfor j = 1:n
    x(j) = pointsEM{j}(1);
    y(j) = pointsEM{j}(2);
end
% xMan = zeros(length(manifold), 1);
% yMan = zeros(length(manifold), 1);
% parfor j = 1:length(manifold)
%     xMan(j) = manifold{j}(1);
%     yMan(j) = manifold{j}(2);
% end
scatter(x, y, 5, colorMap(flagsEM+1,:), 'filled', 'HandleVisibility', 'off')
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
scatter(1-mu, 0, 20, 'w', 'filled', 'DisplayName', "Moon")
plot(cos(linspace(0, 2*pi, 101))*(1-mu), sin(linspace(0, 2*pi, 101)), 'w:', 'DisplayName', "Moon Radius")
% plot(RH/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
axis equal
% ax = axis;
% scatter(xMan, yMan, 5, 'w', 'filled', 'DisplayName', "Stable Manifold")
% axis(ax)
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
title("Periapsis Map $JC_{EM}="+num2str(JC)+"$", 'Interpreter', 'latex')
cb1 = colorbar;  % create the colorbar
colormap(colorMap(1:7,:));
caxis([-0.5 6.5]);
cb1.Ticks = 0:6;
cb1.TickLabels = string(0:6);
cb1.Label.String = 'Periapses';
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'PeriapsisMap_1.png', 'BackgroundColor', 'k')

%% Sun-B1 Periapsis Map
mi = 10;
fig2 = figure("Position", [200 100 1200 750]);
hold on
xSB1 = zeros(n, 1);
ySB1 = zeros(n, 1);
parfor j = 1:n
    xSB1(j) = pointsSB1{mi}{j}(1);
    ySB1(j) = pointsSB1{mi}{j}(2);
end
scatter(xSB1, ySB1, 3, colorMap(flagsSB1{mi}+1,:), 'filled', 'HandleVisibility', 'off')
Earth = plot3DBody("Earth", RE/lstarSB1, [lstar/lstarSB1*cos(moonAngles(mi))*(-mu)+1-muSB1, lstar/lstarSB1*sin(moonAngles(mi))*(-mu), 0]);
set(Earth, 'DisplayName', "Earth")
scatter(lstar/lstarSB1*cos(moonAngles(mi))*(1-mu)+1-muSB1, lstar/lstarSB1*sin(moonAngles(mi))*(1-mu), 20, 'w', 'filled', 'DisplayName', "Moon")
plot(lstar/lstarSB1*cos(linspace(0, 2*pi, 101))*(1-mu)+1-muSB1, lstar/lstarSB1*sin(linspace(0, 2*pi, 101))*(1-mu), 'w:', 'DisplayName', "Moon Radius")
% plot(RH/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
% scatter(xMan.*lstar./lstarSB1+1-muSB1, yMan.*lstar./lstarSB1, 5, 'w', 'filled', 'DisplayName', "Stable Manifold")
axis equal
grid on
xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
title("Periapsis Map $JC_{EM}="+num2str(JC)+"$ $|$ $\theta_{M}="+num2str(moonAngles(mi)*180/pi)+"^\circ$", 'Interpreter', 'latex')
cb2 = colorbar;
colormap(colorMap(1:7,:));
caxis([-0.5 6.5]);
cb2.Ticks = 0:6;
cb2.TickLabels = string(0:6);
cb2.Label.String = 'Periapses';
leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
exportgraphics(fig2, 'PeriapsisMap_2.png', 'BackgroundColor', 'k')

%% Trajectory Plots
% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% B1 = plot3DBody("Earth", RE/lstarSB1, [1-muSB1, 0, 0]);
% set(B1, 'DisplayName', "$B_{1}$")
% plot(lstar/lstarSB1*(cos(linspace(0, 2*pi, 101))*(1-mu))+1-muSB1, lstar/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w:', 'DisplayName', "Moon Radius")
% plot(RH/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstarSB1*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% plot(Traj41.x, Traj41.y, 'c', 'DisplayName', "Trajectory")
% scatter(Traj41.x(1), Traj41.y(1), 50, 'w', '*', 'DisplayName', "S$B_{1}$ Peri.")
% axis equal
% grid on
% xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
% title("Sun-$B_{1}$ Rot.", 'Interpreter', 'latex')
% leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig3, 'PeriapsisMap_3.png', 'BackgroundColor', 'k')
% 
% fig4 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% plot(Traj12.x, Traj12.y, 'c', 'DisplayName', "Trajectory")
% scatter(Traj12.x(1), Traj12.y(1), 50, 'w', '*', 'DisplayName', "S$B_{1}$ Peri.")
% axis equal
% grid on
% xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig4, 'PeriapsisMap_4.png', 'BackgroundColor', 'k')

%% Animation
% xSB1_ani = zeros(length(moonAngles), n);
% ySB1_ani = zeros(length(moonAngles), n);
% parfor j = 1:length(moonAngles)
%     for k = 1:n
%         xSB1_ani(j,k) = pointsSB1{j}{k}(1);
%         ySB1_ani(j,k) = pointsSB1{j}{k}(2);
%     end
% end
% 
% % Video Object and Figure Params
% writerObj = VideoWriter('PeriapsisMap_3_0663_200.mp4', 'MPEG-4');
% writerObj.FrameRate = 3;
% writerObj.Quality = 100;
% open(writerObj);
% fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% fig.Color = 'k';
% grid on
% box on
% 
% % Plot Loop
% for t = 1:length(moonAngles)
%     fig;
% 
%     hold on
%     scatter(xSB1_ani(t,:), ySB1_ani(t,:), 4, colorMap(flagsSB1{t}+1,:), 'filled', 'HandleVisibility', 'off')
%     Earth = plot3DBody("Earth", RE/lstarSB1, [lstar/lstarSB1*cos(moonAngles(t))*(-mu)+1-muSB1, lstar/lstarSB1*sin(moonAngles(t))*(-mu), 0]);
%     set(Earth, 'DisplayName', "Earth")
%     scatter(lstar/lstarSB1*cos(moonAngles(t))*(1-mu)+1-muSB1, lstar/lstarSB1*sin(moonAngles(t))*(1-mu), 20, 'w', 'filled', 'DisplayName', "Moon")
%     plot(lstar/lstarSB1*cos(linspace(0, 2*pi, 101))*(1-mu)+1-muSB1, lstar/lstarSB1*sin(linspace(0, 2*pi, 101))*(1-mu), 'w:', 'DisplayName', "Moon Radius")
%     scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
%     scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
%     axis equal
%     grid on
%     xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
%     ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
%     title("Periapsis Map $JC_{EM}="+num2str(JC)+"$ $|$ $\theta_{M}="+num2str(moonAngles(t)*180/pi)+"^\circ$", 'Interpreter', 'latex')
%     cb2 = colorbar;
%     colormap(colorMap(1:7,:));
%     caxis([-0.5 6.5]);
%     cb2.Ticks = 0:6;
%     cb2.TickLabels = string(0:6);
%     cb2.Label.String = 'Periapses';
%     leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
%     set(gca, 'Color', 'k');
%     view(2)
%     hold off
% 
%     frame = getframe(gcf);
%     writeVideo(writerObj, frame);
%     if t == length(moonAngles)
%         writeVideo(writerObj, frame);
%     end
%     clf
%     grid on
%     box on
% end
% 
% close(writerObj);