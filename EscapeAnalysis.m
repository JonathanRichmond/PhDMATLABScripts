%%% EscapeAnalysis.jl
%%% Jonathan Richmond
%%% C: 22 September 2025
%%% U: 23 September 2025

clear
load('../PhDScripts/Output/EscapeAnalysis.mat')

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

%% Earth-Moon Trajectory Plot
fig1 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
trajs = [Trajectory1, Trajectory2, Trajectory3, Trajectory4, Trajectory5, Trajectory6, Trajectory7];
for j = 1:length(trajs)
    plot3(trajs(j).x(1:end-100), trajs(j).y(1:end-100), trajs(j).z(1:end-100), 'DisplayName', "Escape Traj. "+num2str(j))
end
% plot3(BCR4BPManifold.orbit.x, BCR4BPManifold.orbit.y, BCR4BPManifold.orbit.z, 'b', 'DisplayName', "BCR4BP Orbit")
plot3(CR3BPCompOrbit.x, CR3BPCompOrbit.y, CR3BPCompOrbit.z, 'r', 'DisplayName', "CR3BP Orbit")
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
t2 = tiledlayout(3,1);
t2.TileSpacing = 'compact';
t2.Padding = 'compact';
ax1 = nexttile;
hold on
ps = gobjects(size(d, 1), 1);
for j = 1:size(d, 1)
    ps(j) = plot(d{j}, a_osc_S{j}, 'DisplayName', "$v_{E}="+num2str(v_esc(j))+"$");
end
xline(RH, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
ylim(ax1, [1.3 1.7]*1E8)
grid on
ylabel("$a_{S}$ [km]", 'Interpreter', 'latex')
title("Osculating Elements")
hold off
ax2 = nexttile;
hold on
for j = 1:size(d, 1)
    plot(d{j}, e_osc_S{j}, 'HandleVisibility', 'off')
end
xline(RH, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
ylim(ax2, [0 0.2])
grid on
ylabel("$e_{S}$", 'Interpreter', 'latex')
hold off
ax3 = nexttile;
hold on
for j = 1:size(d, 1)
    plot(d{j}, e_osc_E{j}, 'HandleVisibility', 'off')
end
xline(RH, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
yline(1, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
ylim(ax3, [0 10])
grid on
ylabel("$e_{E}$", 'Interpreter', 'latex')
hold off
xlabel("Earth Distance [km]", 'Interpreter', 'latex')
linkaxes([ax1, ax2, ax3], 'x')
leg2 = legend(ps);
leg2.Layout.Tile = 'east';
leg2.Interpreter = 'latex';
% exportgraphics(fig2, 'EscapeAnalysis_2.png', 'BackgroundColor', 'k')

fig3 = figure("Position", [200 100 1200 750]);
t3 = tiledlayout(3,1);
t3.TileSpacing = 'compact';
t3.Padding = 'compact';
ax1 = nexttile;
hold on
ps = gobjects(size(t, 1), 1);
for j = 1:size(t, 1)
    ps(j) = plot(t{j}, a_osc_S{j}, 'DisplayName', "$v_{E}="+num2str(v_esc(j))+"$");
end
ylim(ax1, [1.3 1.7]*1E8)
grid on
ylabel("$a_{S}$ [km]", 'Interpreter', 'latex')
title("Osculating Elements")
hold off
ax2 = nexttile;
hold on
for j = 1:size(t, 1)
    plot(t{j}, e_osc_S{j}, 'HandleVisibility', 'off')
end
ylim(ax2, [0 0.2])
grid on
ylabel("$e_{S}$", 'Interpreter', 'latex')
hold off
ax3 = nexttile;
hold on
for j = 1:size(t, 1)
    plot(t{j}, e_osc_E{j}, 'HandleVisibility', 'off')
end
yline(1, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
ylim(ax3, [0 10])
grid on
ylabel("$e_{E}$", 'Interpreter', 'latex')
hold off
xlabel("t [EM ndim]", 'Interpreter', 'latex')
linkaxes([ax1, ax2, ax3], 'x')
leg3 = legend(ps);
leg3.Layout.Tile = 'east';
leg3.Interpreter = 'latex';
% exportgraphics(fig3, 'EscapeAnalysis_3.png', 'BackgroundColor', 'k')

fig4 = figure("Position", [200 100 1200 750]);
hold on
for j = 1:size(d, 1)
    plot(a_osc_E{j}, e_osc_E{j}, 'DisplayName', "$v_{E}="+num2str(v_esc(j))+"$")
end
axis([-1E8 1E8 0 2])
grid on
xlabel("$a_{E}$ [km]", 'Interpreter', 'latex')
ylabel("$e_{E}$", 'Interpreter', 'latex')
title("Earth-Moon Rot.", 'Interpreter', 'latex')
leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
hold off
% exportgraphics(fig4, 'EscapeAnalysis_4.png', 'BackgroundColor', 'k')

%% Sun-B1 Trajectory Plot
fig5 = figure("Position", [200 100 1200 750]);
hold on
% Sun = plot3DBody("Sun", RS/lstarSB1, [-muSB1, 0, 0]);
% set(Sun, 'DisplayName', "Sun")
B1 = plot3DBody("Earth", RE/lstarSB1, [1-muSB1, 0, 0]);
set(B1, 'DisplayName', "Earth")
SB1Trajs = [SB1Trajectory1, SB1Trajectory2, SB1Trajectory3, SB1Trajectory4, SB1Trajectory5, SB1Trajectory6, SB1Trajectory7];
for j = 1:length(SB1Trajs)
    plot3(SB1Trajs(j).x(1:end-10), SB1Trajs(j).y(1:end-10), SB1Trajs(j).z(1:end-10), 'DisplayName', "Escape Traj. "+num2str(j))
end
plot3(RH/lstarSB1*cos(linspace(0, 2*pi, 101))+1-muSB1, RH/lstarSB1*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Hills Radius")
axis equal
grid on
xlabel("$x$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
ylabel("$y$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
zlabel("$z$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
title("Sun-$B_{1}$ Rot.", 'Interpreter', 'latex')
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig5, 'EscapeAnalysis_5.png', 'BackgroundColor', 'k')

%% Sun-B1 Hamiltonian
fig6 = figure("Position", [200 100 1200 750]);
t6 = tiledlayout(2,1);
t6.TileSpacing = 'compact';
t6.Padding = 'compact';
ax1 = nexttile;
hold on
colors = ax1.ColorOrder;
nColors = size(colors, 1);
for j = 1:size(t, 1)
    change = [1; find(diff(quad_SB1{j}) ~= 0)+1; length(quad_SB1{j})];
    color = colors(mod(j-1, nColors)+1,:);
    for k = 1:numel(change)-1
        range = change(k):change(k+1);
        quad = quad_SB1{j}(range(1));
        if (quad == 1) || (quad == 3)
            plot(t{j}(range), H_SB1{j}(range), 'Color', color, 'HandleVisibility', 'off')
        else
            plot(t{j}(range), H_SB1{j}(range), 'LineStyle', '--', 'Color', color, 'HandleVisibility', 'off')
        end
    end
end
yline(H_E1, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
grid on
ylabel("$H_{SB_{1}}$ [S$B_{1}$ ndim]", 'Interpreter', 'latex')
title("Sun-$B_{1}$ Hamiltonian", 'Interpreter', 'latex')
hold off
ax2 = nexttile;
hold on
ps = gobjects(size(t, 1), 1);
for j = 1:size(t, 1)
    ps(j) = plot(t{j}, e_osc_E{j}, 'DisplayName', "$v_{E}="+num2str(v_esc(j))+"$");
end
yline(1, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
ylim(ax2, [0 10])
grid on
ylabel("$e_{E}$", 'Interpreter', 'latex')
hold off
xlabel("t [EM ndim]", 'Interpreter', 'latex')
linkaxes([ax1, ax2], 'x')
leg6 = legend(ps);
leg6.Layout.Tile = 'east';
leg6.Interpreter = 'latex';
% exportgraphics(fig6, 'EscapeAnalysis_6.png', 'BackgroundColor', 'k')