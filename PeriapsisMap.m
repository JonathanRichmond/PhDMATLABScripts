%%% PeriapsisMap.jl
%%% Jonathan Richmond
%%% C: 24 October 2025
%%% U: 18 November 2025

clear
% set(0, 'DefaultFigureVisible', 'off');
set(0, 'DefaultFigureVisible', 'on');
load('../PhDScripts/Output//PeriapsisMaps/PeriapsisMap.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0663.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0663_Moon.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0_Moon.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0663_Manifold.mat')
n = length(flagsBCR4BP);

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
colorMap = viridis(6);
colorMap(7,:) = 0.7*[1, 1, 1];
colorMap(8,:) = [1, 0, 0];
colorMap(9,:) = [0, 0, 0];

% peri = 5;
% for j = 1:8
%     if j ~= peri+1
%         colorMap(j,:) = 0.7*[1, 1, 1];
%     end
% end

%% CR3BP Periapsis Map
fig1 = figure("Position", [200 100 1200 750]);
hold on
scatter(xPoints, yPoints, 1.75, colorMap(flagsCR3BP+1,:), 'filled', 'HandleVisibility', 'off')
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
% plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
axis equal
% ax = axis;
% scatter(xManifold, yManifold, 3, 'w', 'filled', 'DisplayName', "Stable Manifold")
% axis(ax)
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
title("Periapsis Map $JC="+num2str(JC)+"$", 'Interpreter', 'latex')
cb1 = colorbar;
colormap(colorMap(1:6,:));
caxis([-0.5 5.5]);
cb1.Ticks = 0:5;
cb1.TickLabels = string(0:5);
cb1.Label.String = 'Periapses';
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'PeriapsisMap_1.png', 'BackgroundColor', 'k')

%% BCR4BP Periapsis Map
t = 1;
fig2 = figure("Position", [200 100 1200 750]);
hold on
scatter(xPoints, yPoints, 1.75, colorMap(flagsBCR4BP{t}+1,:), 'filled', 'HandleVisibility', 'off')
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
% plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
axis equal
% ax = axis;
% scatter(xManifold, yManifold, 3, 'w', 'filled', 'DisplayName', "Stable Manifold")
% axis(ax)
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
title("Periapsis Map $JC="+num2str(JC)+"$ $|$ $\theta_{S}="+num2str(thetaS(t)*180/pi)+"^\circ$", 'Interpreter', 'latex')
cb2 = colorbar;
colormap(colorMap(1:6,:));
caxis([-0.5 5.5]);
cb2.Ticks = 0:5;
cb2.TickLabels = string(0:5);
cb2.Label.String = 'Periapses';
leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig2, 'PeriapsisMap_2.png', 'BackgroundColor', 'k')

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
% plot(RH/lstar*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% plot(PreBurn.x, PreBurn.y, 'Color', 0.7*[1, 1, 1], 'DisplayName', "Old Traj.")
% plot(NoBurn.x, NoBurn.y, 'Color', 0.7*[1, 1, 1], 'HandleVisibility', 'off')
% plot(PostBurn.x, PostBurn.y, 'g', 'DisplayName', "New Traj.")
% scatter(PostBurn.x(1), PostBurn.y(1), 50, 'w', '^', 'filled', 'DisplayName', "$\Delta v="+num2str(deltav*1000, '%.2f')+"$ m/s")
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

%% All Angles Map
% flagsM = cell2mat(flagsBCR4BP(:).');   % each row is flagsBCR4BP{j}
% flagsAllBCR4BP = flagsBCR4BP{1};
% flagsAllBCR4BP(any(flagsM == peri, 2)) = peri;
% xLobe = xPoints(flagsAllBCR4BP == peri);
% yLobe = yPoints(flagsAllBCR4BP == peri);
% area = length(xLobe)/length(flagsAllBCR4BP);
% 
% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xPoints, yPoints, 1.75, colorMap(flagsAllBCR4BP+1,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% % plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
% scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
% axis equal
% grid on
% xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
% title("Periapsis Map (BCR4BP) $JC="+num2str(JC)+"$", 'Interpreter', 'latex')
% cb5 = colorbar;
% colormap(colorMap(1:6,:));
% caxis([-0.5 5.5]);
% cb5.Ticks = 0:5;
% cb5.TickLabels = string(0:5);
% cb5.Label.String = 'Periapses';
% leg5 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig5, 'PeriapsisMap_5.png', 'BackgroundColor', 'k')

%% Animation
% % Video Object and Figure Params
% writerObj = VideoWriter('PeriapsisMap_3_0663_500_Moon.mp4', 'MPEG-4');
% writerObj.FrameRate = 20;
% writerObj.Quality = 100;
% open(writerObj);
% fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% fig.Color = 'k';
% grid on
% box on
% 
% % Plot Loop
% for t = 1:length(thetaS)
%     fig;
% 
%     hold on
%     scatter(xPoints, yPoints, 1.75, colorMap(flagsBCR4BP{t}+1,:), 'filled', 'HandleVisibility', 'off')
%     Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
%     set(Earth, 'DisplayName', "Earth")
%     Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
%     set(Moon, 'DisplayName', "Moon")
%     scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
%     scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
%     axis equal
%     ax = axis;
%     scatter(xManifold, yManifold, 3, 'w', 'filled', 'DisplayName', "Stable Manifold")
%     axis(ax)
%     grid on
%     xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
%     ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
%     title("Periapsis Map $JC="+num2str(JC)+"$ $|$ $\theta_{S}="+num2str(thetaS(t)*180/pi)+"^\circ$", 'Interpreter', 'latex')
%     cb = colorbar;
%     colormap(colorMap(1:6,:));
%     caxis([-0.5 5.5]);
%     cb.Ticks = 0:5;
%     cb.TickLabels = string(0:5);
%     cb.Label.String = 'Periapses';
%     leg = legend('Location', 'bestoutside', 'Interpreter', 'latex');
%     set(gca, 'Color', 'k');
%     view(2)
%     hold off
% 
%     frame = getframe(gcf);
%     writeVideo(writerObj, frame);
%     if t == length(thetaS)
%         writeVideo(writerObj, frame);
%     end
%     clf
%     grid on
%     box on
% end
% 
% close(writerObj);