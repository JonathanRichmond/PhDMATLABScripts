%%% EscapeAnalysis.jl
%%% Jonathan Richmond
%%% C: 22 September 2025
%%% U: 23 September 2025

clear
load('../PhDScripts/Output/EscapeAnalysis.mat')

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
% plot3(BCR4BPManifold.orbit.x, BCR4BPManifold.orbit.y, BCR4BPManifold.orbit.z, 'b', 'DisplayName', "BCR4BP Orbit")
plot3(CR3BPCompOrbit.x, CR3BPCompOrbit.y, CR3BPCompOrbit.z, 'r', 'DisplayName', "CR3BP Orbit")
traj1 = Trajectory2;
plot3(traj1.x, traj1.y, traj1.z, 'b', 'DisplayName', 'Escape Traj. 1')
traj2 = Trajectory3;
plot3(traj2.x, traj2.y, traj2.z, 'r', 'DisplayName', 'Escape Traj. 2')
plot3(RH/lstar*cos(linspace(0, 2*pi, 101))-mu, RH/lstar*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Hills Radius")
axis equal
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
title("Earth-Moon Rot.", 'Interpreter', 'latex')
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'EscapeAnalysis_1.png', 'BackgroundColor', 'k')

%% Osculating Orbital Elements
fig2 = figure("Position", [200 100 1200 750]);
% tiledlayout(2,1)
tiledlayout(3,1)
fig.TileSpacing = 'compact';
fig.Padding = 'compact';
ax1 = nexttile;
hold on
ps = gobjects(size(d, 1), 1);
for j = 1:size(d, 1)
    ps(j) = plot(d{j}, a_osc{j}, 'DisplayName', "$v_{E}="+num2str(v_esc(j))+"$");
end
xline(RH, 'w:', 'LineWidth', 1)
ylim(ax1, [1.3 1.7]*1E8)
grid on
ylabel("$a$ [km]", 'Interpreter', 'latex')
title("Osculating Elements")
hold off
ax2 = nexttile;
hold on
for j = 1:size(d, 1)
    plot(d{j}, e_osc_S{j})
end
xline(RH, 'w:', 'LineWidth', 1)
ylim(ax2, [0 0.2])
grid on
ylabel("$e_{S}$", 'Interpreter', 'latex')
hold off
ax3 = nexttile;
hold on
for j = 1:size(d, 1)
    plot(d{j}, e_osc_E{j})
end
xline(RH, 'w:', 'LineWidth', 1)
yline(1, 'w:')
ylim(ax3, [0 2])
grid on
ylabel("$e_{E}$", 'Interpreter', 'latex')
hold off
xlabel("Earth Distance [km]", 'Interpreter', 'latex')
% linkaxes([ax1, ax2], 'x')
linkaxes([ax1, ax2, ax3], 'x')
leg2 = legend(ps);
leg2.Layout.Tile = 'east';
leg2.Interpreter = 'latex';
% exportgraphics(fig2, 'EscapeAnalysis_2.png', 'BackgroundColor', 'k')