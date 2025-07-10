%%% ManifoldComparison.jl
%%% Jonathan Richmond
%%% C: 10 June 2025
%%% U: 10 July 2025

clear
load('../PhDScripts/Output/ManifoldComparison.mat')

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

RH = 1.496557260502363E6; % Hills radius [km]

%% Trajectory Plot
fig1 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
plot3(BCR4BPManifold.orbit.x, BCR4BPManifold.orbit.y, BCR4BPManifold.orbit.z, 'b', 'DisplayName', "BCR4BP Orbit")
plot3(CR3BPManifold.orbit.x, CR3BPManifold.orbit.y, CR3BPManifold.orbit.z, 'r', 'DisplayName', "CR3BP Orbit")
% for j = 1:BCR4BPManifold.n
%     arc = plotColorLine(BCR4BPManifold.arcs{j}.x, BCR4BPManifold.arcs{j}.y, BCR4BPManifold.arcs{j}.z, BCR4BPManifold.arcs{j}.theta4, 0.005);
%     if j == 1
%         set(arc, 'HandleVisibility', 'on', 'DisplayName', "BCR4BP Manifold")
%     end
% end
% for j = 1:CR3BPManifold.n
%     arc = plot3(CR3BPManifold.arcs{j}.x, CR3BPManifold.arcs{j}.y, CR3BPManifold.arcs{j}.z, 'r', 'LineWidth', 1, 'HandleVisibility', 'off', 'DisplayName', "CR3BP Manifold");
%     if j == 1
%         set(arc, 'HandleVisibility', 'on')
%     end
% end
% for j = 1:pseudoManifold.n
%     V = [pseudoManifold.arcs{j}.x pseudoManifold.arcs{j}.y pseudoManifold.arcs{j}.z];
%     F = [(1:size(V, 1)-1)' (2:size(V, 1))'];
%     arc = patch('Vertices', V, 'Faces', F, 'EdgeColor', 'g', 'EdgeAlpha', 0.5, 'LineWidth', 1, 'FaceColor', 'none', 'HandleVisibility', 'off', 'DisplayName', "Pseudo-Manifolds"); % "P-M ($\theta_{S_{0}}="+num2str(pseudoManifold.theta40*180/pi, '%.0f')+"^{\circ}$)"
%     if j == 1
%         set(arc, 'HandleVisibility', 'on')
%     end
% end
% plot3(RH/lstar*cos(linspace(0, 2*pi, 101))-mu, RH/lstar*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Hills Radius")
axis equal
axis
% axisLimits = [-0.22 1.52 -0.57 0.57];
grid on
% ZVC = getCR3BPZVC(mu, CR3BPManifold.orbit.JC, axisLimits, 250, 'w--');
% set(ZVC, 'DisplayName', "CR3BP ZVC")
% axis(axisLimits)
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
title("Earth-Moon Rot.", 'Interpreter', 'latex')
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
phasemap
% phasebar('deg', 'Location', 'northwest', 'Size', 0.275);
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'ManifoldComparison_1.png', 'BackgroundColor', 'k')

%% Metric Plot
fig2 = figure("Position", [200 100 1200 750]);
hold on
scatter(HillsTOFBCR4BP, vEscBCR4BP, 20, [1 0 1], 'filled', 'DisplayName', "BCR4BP Manifold")
scatter(HillsTOFCR3BP, vEscCR3BP, 20, 'r', 'filled', 'DisplayName', "CR3BP Manifold")
scatter(HillsTOFPseudo, vEscPseudo, 20, 'g', 'filled', 'DisplayName', "Pseudo-Manifolds")
axis([0 9.99*pi 0 1.2])
grid on
xlabel("$TOF$ [EM ndim]", 'Interpreter', 'latex')
ylabel("ECI $v_{esc}$ [EM ndim]", 'Interpreter', 'latex')
title("Hills Escape Velocity", 'Interpreter', 'latex')
legend('Location', 'bestoutside', 'Interpreter', 'latex');
hold off
% exportgraphics(fig2, 'ManifoldComparison_2.png', 'BackgroundColor', 'k')