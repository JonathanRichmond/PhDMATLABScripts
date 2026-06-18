%%% EscapeCR3BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 16 June 2026
%%% U: 18 June 2026

clear

%% Import Map Data
volumeData = load('../PhDScripts/Output/ApseMaps/CR3BPJCVolume_1_peri_pro_500_2.8_3.1.mat');

mapsData = load('../PhDScripts/Output/ApseMaps/TestMap.mat');
map = mapsData.map;
primary = map.primary;
switch map.grade
    case "pro"
        grade = "prograde";
    case "retro"
        grade = "retrograde";
end
switch map.apse
    case "peri"
        apse = "periapsis";
    case "apo"
        apse = "apoapsis";
end
JC = map.JC;
disp(primary+"-centered "+grade+" "+apse+" map: JC = "+JC)

n = size(map.flags, 1);
xGrid = map.q(1,:);
yGrid = map.q(2,:);
flags = reshape(map.flags, 1, n^2);

%% Earth-Moon Data
gmE = 3.9860043543609593E5; % Earth gravitational parameter [km^3/s^2]
mE = gmE/6.67384E-20; % Earth mass [kg]
RE = 6.371008366666666E3; % Earth radius [km]

gmm = 4.9028000661637961E3; % Moon gravitational parameter [km^3/s^2]
mm = gmm/6.67384E-20; % Moon mass [kg]
Rm = 1.7374000000000003E3; % Moon radius [km]

muEM = gmm/(gmE+gmm); % Mass ratio
mstarEM = (gmE+gmm)/6.67384E-20; % Characteristic mass [kg]
lstarEM = 3.8474799201129237E5; % Characteristic length [km]
tstarEM = sqrt(lstarEM^3/(gmE+gmm)); % Characteristic time [s]

g1EM = muEM; % Initial guess
delg1 = 1;
while abs(delg1) > eps
    f = ((1-muEM)/((1-g1EM)^2))-(muEM/(g1EM^2))-1+muEM+g1EM;
    fprime = ((2*(1-muEM))/((1-g1EM)^3))+((2*muEM)/(g1EM^3))+1;
    g1EMnew = g1EM-(f/fprime);
    delg1 = g1EMnew-g1EM;
    g1EM = g1EMnew;
    i = i+1;
end
a1EM = 1-muEM-g1EM;

g2EM = muEM; % Initial guess
delg2 = 1;
while abs(delg2) > eps
    f = ((1-muEM)/((1+g2EM)^2))+(muEM/(g2EM^2))-1+muEM-g2EM;
    fprime = ((-2*(1-muEM))/((1+g2EM)^3))-((2*muEM)/(g2EM^3))-1;
    g2EMnew = g2EM-(f/fprime);
    delg2 = g2EMnew-g2EM;
    g2EM = g2EMnew;
    i = i+1;
end
a2EM = 1-muEM+g2EM;

g3EM = muEM; %Initial guess
delg3 = 1;
while abs(delg3) > eps
    f = ((1-muEM)/(g3EM^2))+(muEM/((1+g3EM)^2))-muEM-g3EM;
    fprime = ((-2*(1-muEM))/(g3EM^3))-((2*muEM)/((1+g3EM)^3))-1;
    g3EMnew = g3EM-(f/fprime);
    delg3 = g3EMnew-g3EM;
    g3EM = g3EMnew;
    i = i+1;
end
a3EM = -1*muEM-g3EM;

a45EM = 0.5-muEM;
b4EM = sqrt(3)/2;
b5EM = -b4EM;

%% Sun-Earth Data
gmS = 1.327124400419393E11; % Sun gravitational parameter [km^3/s^2]
mS = gmS/6.67384E-20; % Sun mass [kg]
RS = 6.96E5; % Sun radius [km]

muSE = gmE/(gmS+gmE); % Mass ratio
mstarSE = (gmS+gmE)/6.67384E-20; % Characteristic mass [kg]
lstarSE = 1.4959789217545033E8; % Characteristic length [km]
tstarSE = sqrt(lstarSE^3/(gmS+gmE)); % Characteristic time [s]

g1SE = muSE; % Initial guess
delg1 = 1;
while abs(delg1) > eps
    f = ((1-muSE)/((1-g1SE)^2))-(muSE/(g1SE^2))-1+muSE+g1SE;
    fprime = ((2*(1-muSE))/((1-g1SE)^3))+((2*muSE)/(g1SE^3))+1;
    g1SEnew = g1SE-(f/fprime);
    delg1 = g1SEnew-g1SE;
    g1SE = g1SEnew;
    i = i+1;
end
a1SE = 1-muSE-g1SE;

g2SE = muSE; % Initial guess
delg2 = 1;
while abs(delg2) > eps
    f = ((1-muSE)/((1+g2SE)^2))+(muSE/(g2SE^2))-1+muSE-g2SE;
    fprime = ((-2*(1-muSE))/((1+g2SE)^3))-((2*muSE)/(g2SE^3))-1;
    g2SEnew = g2SE-(f/fprime);
    delg2 = g2SEnew-g2SE;
    g2SE = g2SEnew;
    i = i+1;
end
a2SE = 1-muSE+g2SE;

g3SE = muSE; %Initial guess
delg3 = 1;
while abs(delg3) > eps
    f = ((1-muSE)/(g3SE^2))+(muSE/((1+g3SE)^2))-muSE-g3SE;
    fprime = ((-2*(1-muSE))/(g3SE^3))-((2*muSE)/((1+g3SE)^3))-1;
    g3SEnew = g3SE-(f/fprime);
    delg3 = g3SEnew-g3SE;
    g3SE = g3SEnew;
    i = i+1;
end
a3SE = -1*muSE-g3SE;

a45SE = 0.5-muSE;
b4SE = sqrt(3)/2;
b5SE = -b4SE;

RSoIE = 0.09877*lstarSE; % Earth sphere of influence radius [km]

%% Colormap
colorMap = viridis(6); % Escape
colorMap(7,:) = 0.3*[1, 1, 1]; % Capture
colorMap(8,:) = [1, 0, 0]; % Impact
% colorMap(8,:) = 0.3*[1, 1, 1]; % Impact
colorMap(9,:) = [0, 0, 0]; % Invalid apse
colorMap(10,:) = [1, 1, 1]; % ZVC

%% Map
fig1 = figure("Position", [200 100 1200 750]);
hold on
scatter(xGrid, yGrid, 1.75, colorMap(flags+1,:), 'filled', 'HandleVisibility', 'off')
Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
set(Moon, 'DisplayName', "Moon")
scatter(nan, nan, 20, colorMap(7,:), 'filled', 'DisplayName', "Capture")
scatter(nan, nan, 20, colorMap(8,:), 'filled', 'DisplayName', "Impact")
scatter(nan, nan, 20, colorMap(10,:), 'filled', 'DisplayName', "ZVC")
axis equal
axis([-1.25 1.25 -1.25 1.25])
% axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
colormap(colorMap(1:6,:))
cb1 = colorbar;
clim([-0.5 5.5])
cb1.Ticks = 0:5;
cb1.TickLabels = [string(0:4), "5+"];
cb1.Label.String = 'Periapses';
ylabel(cb1, "Periapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
cb1.Label.Position = cb1.Label.Position+[-1.5 3.1 0];
leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
drawnow;
set(leg1.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
set(gca, 'Color', 'k');
view(2)
hold off
ax1 = gca;
ax1.SortMethod = 'childorder';
% exportgraphics(fig1, 'EscapeCR3BP_1.png', 'BackgroundColor', 'k')

%% Propagators
odeCR3BPEM = @(t,r) ODE_CR3BP(t, r, muEM);
odeOpts = odeset('RelTol', 1E-12, 'AbsTol', 1E-12);

%% Test Trajectory
% xSample = 0.308116;
% ySample = -0.388277;
% idx = find((abs(xGrid-xSample) < 1E-5) & (abs(yGrid-ySample) < 1E-5));
% q = map.q(:,idx);
% disp("Sample IC: ["+q(1)+", "+q(2)+", "+q(3)+", "+q(4)+", "+q(5)+", "+q(6)+"]")
% tau = 0.2*pi;
% sol = ode89(odeCR3BPEM, [0 tau], q, odeOpts);
% 
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% p21 = plot3WithArrows(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', colorMap(flags(idx)+1,:));
% set(p21, 'DisplayName', "Sample Traj.")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% % axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg2.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax2 = gca;
% ax2.SortMethod = 'childorder';
% % exportgraphics(fig2, 'EscapeCR3BP_2.png','BackgroundColor', 'k')