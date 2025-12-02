%%% PeriapsisMap.jl
%%% Jonathan Richmond
%%% C: 24 October 2025
%%% U: 2 December 2025

clear
% set(0, 'DefaultFigureVisible', 'off');
set(0, 'DefaultFigureVisible', 'on');
% load('../PhDScripts/Output//PeriapsisMaps/PeriapsisMap.mat')
load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0663.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0663_Moon.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0_Moon.mat')
% load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMap_3_0663_Manifold.mat')
load('../PhDScripts/Output/PeriapsisMaps/PeriapsisMaps_JC.mat')
n = length(flagsBCR4BP);
nJC = size(flagsCR3BPJC, 1);

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
% colorMap(8,:) = [1, 0, 0];
colorMap(8,:) = 0.7*[1, 1, 1];
colorMap(9,:) = [0, 0, 0];

% peri = 0;
% for j = 1:8
%     if j ~= peri+1
%         colorMap(j,:) = 0.7*[1, 1, 1];
%     end
% end

%% CR3BP Periapsis Map
% fig1 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xPoints, yPoints, 1.75, colorMap(flagsCR3BP+1,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% % plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% % scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
% % scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
% axis([-1.25 1.25 -1.25 1.25])
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% % ax = axis;
% % scatter(xManifold, yManifold, 3, 'w', 'filled', 'DisplayName', "Stable Manifold")
% % axis(ax)
% xlabel("x [EM ndim]")
% ylabel("y [EM ndim]")
% title("Periapsis Map JC = "+num2str(JC))
% cb1 = colorbar;
% colormap(colorMap(1:6,:));
% clim([-0.5 5.5]);
% cb1.Ticks = 0:5;
% cb1.TickLabels = string(0:5);
% cb1.Label.String = 'Periapses';
% % leg1 = legend('Location', 'bestoutside');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% exportgraphics(fig1, 'PeriapsisMap_1.png', 'BackgroundColor', 'k')

%% BCR4BP Periapsis Map
% t = 1;
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xPoints, yPoints, 1.75, colorMap(flagsBCR4BP{t}+1,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% % plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% % scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
% % scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
% axis([-1.25 1.25 -1.25 1.25])
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% % ax = axis;
% % scatter(xManifold, yManifold, 3, 'w', 'filled', 'DisplayName', "Stable Manifold")
% % axis(ax)
% xax = xlabel("x [EM ndim]");
% ylabel("y [EM ndim]")
% title("Periapsis Map JC = "+num2str(JC)+" | \theta_{S} = "+num2str(thetaS(t)*180/pi)+"^\circ")
% cb2 = colorbar;
% colormap(colorMap(1:6,:));
% clim([-0.5 5.5]);
% cb2.Ticks = 0:5;
% cb2.TickLabels = string(0:5);
% cb2.Label.String = 'Periapses';
% % leg2 = legend('Location', 'bestoutside');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig2, 'PeriapsisMap_2.png', 'BackgroundColor', 'k')

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
% plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
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

%% Persistent BCR4BP Escape Map
% flagsM = cell2mat(flagsBCR4BP(:).');
% % percentPeri = sum(flagsM == peri, 2)*100/size(flagsM, 2);
% percentPeri = sum(flagsM < 6, 2)*100/size(flagsM, 2);
% flagsAllBCR4BP = flagsBCR4BP{1};
% flagsAllBCR4BP(flagsAllBCR4BP < 6) = 6;
% % flagsAllBCR4BP(percentPeri >= 33) = peri;
% % xLobe = xPoints(flagsAllBCR4BP == peri);
% % yLobe = yPoints(flagsAllBCR4BP == peri);
% flagsAllBCR4BP(percentPeri >= 33) = 0;
% xLobe = xPoints(flagsAllBCR4BP == 0);
% yLobe = yPoints(flagsAllBCR4BP == 0);
% areaBCR4BP = length(xLobe)*100/length(flagsAllBCR4BP);

% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xPoints, yPoints, 1.75, colorMap(flagsAllBCR4BP+1,:), 'filled', 'HandleVisibility', 'off')
% scatter(xLobe, yLobe, 1.75, 'b', 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% % plot(RH/lstar*cos(linspace(0, 2*pi, 101)), RH/lstar*sin(linspace(0, 2*pi, 101)), 'w--', 'DisplayName', "Hills Radius")
% % scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
% % scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
% axis([-1.25 1.25 -1.25 1.25])
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% xlabel("x [EM ndim]")
% ylabel("y [EM ndim]")
% title("Periapsis Map (BCR4BP) JC = "+num2str(JC))
% % cb5 = colorbar;
% % colormap(colorMap(1:6,:));
% % clim([-0.5 5.5]);
% % cb5.Ticks = 0:5;
% % cb5.TickLabels = string(0:5);
% % cb5.Label.String = 'Periapses';
% % leg5 = legend('Location', 'bestoutside');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig5, 'PeriapsisMap_5.png', 'BackgroundColor', 'k')

%% JC Volume
% flagsGrid = reshape(flagsCR3BPJC.', sqrt(size(flagsCR3BPJC, 2)), sqrt(size(flagsCR3BPJC, 2)), size(flagsCR3BPJC, 1));
% flagsGrid = permute(flagsGrid, [2 1 3]);
% V = uint8(flagsGrid(1:end, 1:end, 1:end)+1);
% alphas = zeros(size(colorMap, 1), 1);
% alphas(1:6) = 1;
% v = volshow(V, 'Colormap', colorMap, 'Alphamap', alphas);
% p = v.Parent;
% % Do manually
% % fig6 = getframe(p);
% % imwrite(fig6.cdata, 'PeriapsisMap_6.png');

%% Area Plot
% areaCR3BP = zeros(1, 6);
% areaBCR4BP = zeros(n, 6);
% for k = 1:6
%     if k == 1
%         areaCR3BP(k) = sum(flagsCR3BP == k-1)*100/length(flagsCR3BP)-33;
%     else
%         areaCR3BP(k) = sum(flagsCR3BP == k-1)*100/length(flagsCR3BP);
%     end
%     for j = 1:n
%         if k == 1
%             areaBCR4BP(j,k) = sum(flagsBCR4BP{j} == k-1)*100/length(flagsBCR4BP{j})-33;
%         else
%             areaBCR4BP(j,k) = sum(flagsBCR4BP{j} == k-1)*100/length(flagsBCR4BP{j});
%         end
%     end
% end
% 
% fig7 = figure("Position", [200 100 1200 750]);
% hold on
% for k = 1:6
%     hy = yline(areaCR3BP(k), '--', 'Color', colorMap(k,:), 'LineWidth', 2, 'DisplayName', "CR3BP", 'HandleVisibility', 'off');
%     p7 = plot(0:n-1, areaBCR4BP(:,k)', 'Color', colorMap(k,:), 'DisplayName', "BCR4BP", 'HandleVisibility', 'off');
%     if k == 1
%         set(hy, 'HandleVisibility', 'on');
%         set(p7, 'HandleVisibility', 'on');
%     end
% end
% ax = axis;
% axis([0 360 0 2])
% grid on
% xlabel("\theta_{S} [deg]")
% ylabel("Total Area [%]")
% title("JC = "+num2str(JC))
% leg7 = legend('Location', 'bestoutside');
% set(gca, 'Color', 'k');
% hold off
% % exportgraphics(fig7, 'PeriapsisMap_7.png', 'BackgroundColor', 'k')

%% Delta v Map
% % target = [-0.398297 -0.207916];
% % sample = find(vecnorm([xPoints(:) yPoints(:)] - target, 2, 2) < 1e-5, 1);
% JCidx = find(JCs-JC < 0, 1);
% newJCs = 1.995*ones(length(xPoints), 1);
% mask = flagsCR3BPJC(JCidx:end,:) < 6;
% [ridx, cidx] = find(mask);
% [~, ri] = unique(cidx, 'first');
% newJCs(cidx(ri)) = JCs(JCidx - 1 + ridx(ri));
% r13 = sqrt((xPoints+mu).^2+yPoints.^2);
% r23 = sqrt((xPoints-1+mu).^2+yPoints.^2);
% Omegas = (1-mu)./r13+mu./r23+0.5.*(xPoints.^2+yPoints.^2);
% DeltavMap = zeros(length(xPoints), 1);
% parfor j = 1:length(xPoints)
%     if flagsAllBCR4BP(j) == 0
%         continue
%     end
%     newJC = newJCs(j);
%     DeltavMap(j) = (sqrt(2*Omegas(j)-newJC)-sqrt(2*Omegas(j)-JC))*lstar*1000/tstar; % m/s
% end
% cplxDeltav = abs(imag(DeltavMap)) > 1E-12;
% realDeltav = ~cplxDeltav;
% 
% fig8 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xPoints(realDeltav), yPoints(realDeltav), 1.75, DeltavMap(realDeltav), 'filled', 'HandleVisibility', 'off')
% scatter(xPoints(cplxDeltav), yPoints(cplxDeltav), 1.75, 'w', 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% axis([-1.25 1.25 -1.25 1.25])
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% xlabel("x [EM ndim]")
% ylabel("y [EM ndim]")
% title("Periapsis Escape \Deltav Map JC = "+num2str(JC))
% cb8 = colorbar;
% colormap(turbo);
% clim([min(DeltavMap(realDeltav)) max(DeltavMap(realDeltav))]);
% lbls8 = string(0:100:900);
% lbls8(end) = ">900";
% cb8.TickLabels = lbls8;
% cb8.Label.String = '\Deltav [m/s]';
% % leg8 = legend('Location', 'bestoutside');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig8, 'PeriapsisMap_8.png', 'BackgroundColor', 'k')

%% Sun Angle Animation
% % Video Object and Figure Params
% writerObj = VideoWriter('PeriapsisMap_3_0663_500_Moon.mp4', 'MPEG-4');
% writerObj.FrameRate = 20;
% writerObj.Quality = 100;
% open(writerObj);
% fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% fig.Color = 'k';
% box on
% 
% % Plot Loop
% for t = 1:length(thetaS)
%     hold on
%     scatter(xPoints, yPoints, 1.75, colorMap(flagsBCR4BP{t}+1,:), 'filled', 'HandleVisibility', 'off')
%     % Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
%     % set(Earth, 'DisplayName', "Earth")
%     Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
%     % set(Moon, 'DisplayName', "Moon")
%     % scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
%     % scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
%     scale = 0.03*[-1 1 -1 1];
%     % axis(scale+[-1.25 1.25 -1.25 1.25])
%     axis(scale+[1-mu-0.3 1-mu+0.3 -0.3 0.3])
%     axis equal
%     % axis(scale+[-1.25 1.25 -1.25 1.25])
%     axis(scale+[1-mu-0.3 1-mu+0.3 -0.3 0.3])
%     xlabel("x [EM ndim]")
%     ylabel("y [EM ndim]")
%     title("Periapsis Map JC = "+num2str(JC)+" | \theta_{S} = "+num2str(thetaS(t)*180/pi)+"^\circ")
%     cb = colorbar;
%     colormap(colorMap(1:6,:));
%     clim([-0.5 5.5]);
%     cb.Ticks = 0:5;
%     cb.TickLabels = string(0:5);
%     cb.Label.String = 'Periapses';
%     leg = legend('Location', 'bestoutside');
%     leg.Visible = 'off';
%     set(gca, 'Color', 'k');
%     view(2)
%     hold off
% 
%     drawnow;
%     frame = getframe(gcf);
%     writeVideo(writerObj, frame);
%     if t == length(thetaS)
%         writeVideo(writerObj, frame);
%     end
%     clf
%     box on
% end
% 
% close(writerObj);

%% JC Animation
% % Video Object and Figure Params
% writerObj = VideoWriter('PeriapsisMap_JC.mp4', 'MPEG-4');
% writerObj.FrameRate = 20;
% writerObj.Quality = 100;
% open(writerObj);
% fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% fig.Color = 'k';
% box on
% 
% % Plot Loop
% for t = 1:size(flagsCR3BPJC, 1)
%     hold on
%     scatter(xPoints, yPoints, 1.75, colorMap(flagsCR3BPJC(t,:)+1,:), 'filled', 'HandleVisibility', 'off')
%     Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
%     % set(Earth, 'DisplayName', "Earth")
%     Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
%     % set(Moon, 'DisplayName', "Moon")
%     % scatter(nan, nan, 20, 0.7*[1, 1, 1], 'filled', 'DisplayName', "Capture")
%     % scatter(nan, nan, 20, [1, 0, 0], 'filled', 'DisplayName', "Impact")
%     scale = 1.1*[-1 1 -1 1];
%     axis(scale+[-1.25 1.25 -1.25 1.25])
%     axis equal
%     axis(scale+[-1.25 1.25 -1.25 1.25])
%     xlabel("x [EM ndim]")
%     ylabel("y [EM ndim]")
%     title("Periapsis Map JC = "+num2str(JCs(t), '%.4f'))
%     cb = colorbar;
%     colormap(colorMap(1:6,:));
%     clim([-0.5 5.5]);
%     cb.Ticks = 0:5;
%     cb.TickLabels = string(0:5);
%     cb.Label.String = 'Periapses';
%     leg = legend('Location', 'bestoutside');
%     leg.Visible = 'off';
%     set(gca, 'Color', 'k');
%     view(2)
%     hold off
% 
%     drawnow;
%     frame = getframe(gcf);
%     writeVideo(writerObj, frame);
%     clf
%     box on
% end
% 
% close(writerObj);