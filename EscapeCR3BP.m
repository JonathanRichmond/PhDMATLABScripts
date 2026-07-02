%%% EscapeCR3BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 16 June 2026
%%% U: 2 July 2026

clear

%% Import Map Data
mapsData = load('../PhDScripts/Output/ApseMaps/CR3BP_1_peri_pro_500_3.0663.mat');
mapFields = fieldnames(mapsData);
map = mapsData.(mapFields{1});
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

%% Sun-Mars Data
gmM = 4.282837362069909E4; % Mars gravitational parameter [km^3/s^2]
mM = gmM/6.67384E-20; % Mars mass [kg]
RM = 3.3895266666666666E3; % Mars radius [km]

muSM = gmM/(gmS+gmM); % Mass ratio
mstarSM = (gmS+gmM)/6.67384E-20; % Characteristic mass [kg]
lstarSM = 2.2794082723873448E8; % Characteristic length [km]
tstarSM = sqrt(lstarSM^3/(gmS+gmM)); % Characteristic time [s]

g1SM = muSM; % Initial guess
delg1 = 1;
while abs(delg1) > eps
    f = ((1-muSM)/((1-g1SM)^2))-(muSM/(g1SM^2))-1+muSM+g1SM;
    fprime = ((2*(1-muSM))/((1-g1SM)^3))+((2*muSM)/(g1SM^3))+1;
    g1SMnew = g1SM-(f/fprime);
    delg1 = g1SMnew-g1SM;
    g1SM = g1SMnew;
    i = i+1;
end
a1SM = 1-muSM-g1SM;

g2SM = muSM; % Initial guess
delg2 = 1;
while abs(delg2) > eps
    f = ((1-muSM)/((1+g2SM)^2))+(muSM/(g2SM^2))-1+muSM-g2SM;
    fprime = ((-2*(1-muSM))/((1+g2SM)^3))-((2*muSM)/(g2SM^3))-1;
    g2SMnew = g2SM-(f/fprime);
    delg2 = g2SMnew-g2SM;
    g2SM = g2SMnew;
    i = i+1;
end
a2SM = 1-muSM+g2SM;

g3SM = muSM; %Initial guess
delg3 = 1;
while abs(delg3) > eps
    f = ((1-muSM)/(g3SM^2))+(muSM/((1+g3SM)^2))-muSM-g3SM;
    fprime = ((-2*(1-muSM))/(g3SM^3))-((2*muSM)/((1+g3SM)^3))-1;
    g3SMnew = g3SM-(f/fprime);
    delg3 = g3SMnew-g3SM;
    g3SM = g3SMnew;
    i = i+1;
end
a3SM = -1*muSM-g3SM;

a45SM = 0.5-muSM;
b4SM = sqrt(3)/2;
b5SM = -b4SM;

%% Colormap
colorMap = viridis(6); % Escape
colorMap(7,:) = [0.78, 0.72, 0.66]; % Capture
colorMap(8,:) = [1, 0, 0]; % Impact
% colorMap(8,:) = [0.78, 0.72, 0.66]; % Impact
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
ylabel(cb1, "Periapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% ylabel(cb1, "Apoapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
cb1.Label.Position = cb1.Label.Position+[-2.3 3.1 0];
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
% xSample = -0.0876754;
% ySample = -0.258016;
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
% scatter3(solOrig.y(1,1), solOrig.y(2,1), solOrig.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
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

%% Choose Sample Trajectory
xSample = -0.0826653;
ySample = -0.222946;
idx = find((abs(xGrid-xSample) < 1E-5) & (abs(yGrid-ySample) < 1E-5))

%% Import Escape Analysis Data
analysisData = load('../PhDScripts/Output/EscapeAnalysisCR3BP.mat');
esc0q0s = analysisData.esc0q0;
esc0tfs = analysisData.esc0tf;
esc0Es = analysisData.esc0E;
esc0n = length(esc0Es);
esc1q0s = analysisData.esc1q0;
esc1tfs = analysisData.esc1tf;
esc1Es = analysisData.esc1E;
esc1n = length(esc1Es);
esc2q0s = analysisData.esc2q0;
esc2tfs = analysisData.esc2tf;
esc2Es = analysisData.esc2E;
esc2n = length(esc2Es);

Deltav1s = analysisData.Deltav1s;
escE1s = analysisData.EscapeEs;

norms0 = vecnorm(esc0q0s(1:2,:));
mask0 = norms0 < 1;
Es_filt0 = esc0Es(mask0);
q0s_filt0 = esc0q0s(:, mask0);
tfs_filt0 = esc0tfs(mask0);
nFilt0 = sum(mask0);
clear esc0q0s esc0tfs esc0Es mask0
norms1 = vecnorm(esc1q0s(1:2,:));
mask1 = norms1 < 1;
Es_filt1 = esc1Es(mask1);
q0s_filt1 = esc1q0s(:, mask1);
tfs_filt1 = esc1tfs(mask1);
nFilt1 = sum(mask1);
clear esc1q0s esc1tfs esc1Es mask1
norms2 = vecnorm(esc2q0s(1:2,:));
mask2 = norms2 < 1;
Es_filt2 = esc2Es(mask2);
q0s_filt2 = esc2q0s(:, mask2);
tfs_filt2 = esc2tfs(mask2);
nFilt2 = sum(mask2);
clear esc2q0s esc2tfs esc2Es mask2
E_min = min([Es_filt0; Es_filt1; Es_filt2]);
E_max = max([Es_filt0; Es_filt1; Es_filt2]);

%% Escape Analysis
DeltavHs = nan(length(escE1s), 1);
r2 = a1SM*lstarSM;
parfor j = 1:length(escE1s)
    E = escE1s(j);
    DeltavHs(j) = abs(sqrt((8*E^2*r2)/(gmS-2*E*r2))-sqrt(-2*E))+abs(sqrt(gmS/r2)-sqrt((2*gmS^2)/(r2*(gmS-2*E*r2))));
end

fig8 = figure("Position", [200 100 1200 750]);
hold on
% scatter(Deltav1s.*1000.*lstarEM./tstarEM, escE1s, 20, 'filled', 'HandleVisibility', 'off')
% scatter(Deltav1s.*1000.*lstarEM./tstarEM, DeltavHs, 20, 'filled', 'HandleVisibility', 'off')
scatter(Deltav1s.*1000.*lstarEM./tstarEM, Deltav1s.*lstarEM./tstarEM+DeltavHs, 20, 'filled', 'HandleVisibility', 'off')
xlabel("$\Delta v_{1}$ [m/s]", 'Interpreter', 'latex')
% ylabel("$\mathcal{E}_{esc}$ [km$^{2}$/s$^{2}$]", 'Interpreter', 'latex')
% ylabel("$\Delta v_{H}$ [km/s]", 'Interpreter', 'latex')
ylabel("Total $\Delta v$ [km/s]", 'Interpreter', 'latex')
title("Maneuver Optimization", 'Interpreter', 'latex')
set(gca, 'Color', 'k');
view(2)
hold off
ax8 = gca;
ax8.SortMethod = 'childorder';
% exportgraphics(fig8, 'EscapeCR3BP_8.png', 'BackgroundColor', 'k')

qOrig = map.q(:,idx);
disp("Original IC: ["+qOrig(1)+", "+qOrig(2)+", "+qOrig(3)+", "+qOrig(4)+", "+qOrig(5)+", "+qOrig(6)+"]")
tauOrig = 3*pi;
solOrig = ode89(odeCR3BPEM, [0 tauOrig], qOrig, odeOpts);
[~, optIdx] = min(DeltavHs);
vOrig = norm(qOrig(4:5));
vhat = qOrig(4:5)./vOrig;
qAssist = qOrig;
Deltav1 = Deltav1s(optIdx);
% Deltav1 = 239.024*tstarEM/1000/lstarEM;
qAssist(4:5) = (vOrig+Deltav1).*vhat;
disp("Delta-v: "+Deltav1*1000*lstarEM/tstarEM+" m/s")
disp("Assisted IC: ["+qAssist(1)+", "+qAssist(2)+", "+qAssist(3)+", "+qAssist(4)+", "+qAssist(5)+", "+qAssist(6)+"]")
tauAssist = 12*pi;
solAssist = ode89(odeCR3BPEM, [0 tauAssist], qAssist, odeOpts);

fig9 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
set(Moon, 'DisplayName', "Moon")
scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
scatter3(solOrig.y(1,1), solOrig.y(2,1), solOrig.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
p91 = plot3WithArrows(solOrig.y(1,:), solOrig.y(2,:), solOrig.y(3,:), 'Color', colorMap(flags(idx)+1,:));
set(p91, 'DisplayName', "Original Traj.")
p92 = plot3WithArrows(solAssist.y(1,:), solAssist.y(2,:), solAssist.y(3,:), 'g');
set(p92, 'DisplayName', "Assisted Traj.")
axis equal
% axis([-1.25 1.25 -1.25 1.25])
axis([-3 3 -3 3])
grid on
xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
title("Earth-Moon Rot.", 'Interpreter', 'latex')
leg9 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
drawnow;
set(leg9.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
set(gca, 'Color', 'k');
view(2)
hold off
ax9 = gca;
ax9.SortMethod = 'childorder';
% exportgraphics(fig9, 'EscapeCR3BP_9.png','BackgroundColor', 'k')

%% Escape Analysis Figure
% colors = nebula(1000);
% pointColors0 = zeros(nFilt0, 3);
% parfor j = 1:nFilt0
%     pointColors0(j,:) = getColor(colors, Es_filt0(j), [E_min, E_max]);
% end
% pointColors1 = zeros(nFilt1, 3);
% parfor j = 1:nFilt1
%     pointColors1(j,:) = getColor(colors, Es_filt1(j), [E_min, E_max]);
% end
% pointColors2 = zeros(nFilt2, 3);
% parfor j = 1:nFilt2
%     pointColors2(j,:) = getColor(colors, Es_filt2(j), [E_min, E_max]);
% end
% 
% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% color = colorMap(1,:);
% for j = 1:nFilt2
%     q = q0s_filt2(:,j);
%     tau = tfs_filt2(j);
%     [~, rout] = mexCR3BP(q, [0 tau], muEM, 1E-12, 1E-12, 1E-10);
%     % plot3(rout(:,1), rout(:,2), rout(:,3), 'Color', [0.5, 0.5, 0.5, 0.1], 'LineWidth', 1, 'HandleVisibility', 'off')
%     plot3(rout(:,1), rout(:,2), rout(:,3), 'Color', [pointColors2(j,:), 0.1], 'LineWidth', 1, 'HandleVisibility', 'off')   
% end
% % scatter3(q0s_filt0(1,:), q0s_filt0(2,:), q0s_filt0(3,:), 1.75, color, 'filled', 'HandleVisibility', 'off')
% % scatter3(q0s_filt0(1,:), q0s_filt0(2,:), q0s_filt0(3,:), 1.75, pointColors0, 'filled', 'HandleVisibility', 'off')
% % scatter3(q0s_filt1(1,:), q0s_filt1(2,:), q0s_filt1(3,:), 1.75, pointColors1, 'filled', 'HandleVisibility', 'off')
% scatter3(q0s_filt2(1,:), q0s_filt2(2,:), q0s_filt2(3,:), 1.75, pointColors2, 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% colormap(nebula)
% cb5 = colorbar;
% clim([E_min E_max])
% ylabel(cb5, "$\mathcal{E}_{esc}$ [km$^{2}$/s^${2}$]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb5.Label.Position = cb5.Label.Position+[-3 5.5 0];
% leg5 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% leg5.Position = leg5.Position+[0.08 0 0 0];
% drawnow;
% set(leg5.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax5 = gca;
% ax5.SortMethod = 'childorder';
% % exportgraphics(fig5, 'EscapeCR3BP_5.png','BackgroundColor', 'k')

%% Import JC Volume Data
volumeFile = 'E:/ApseMapData/CR3BPJCVolume_1_peri_pro_500_0.9_3.15.mat';
volumeDataFile = 'CR3BPJCVolume_1_peri_pro_500_0.9_3.15.mat';

% volumeFields = who('-file', volumeFile);
% nVolume = length(volumeFields);
% sampleData = load(volumeFile, volumeFields{1});
% sampleFields = fieldnames(sampleData);
% sampleMap = sampleData.(sampleFields{1});
% nSample = size(sampleMap.flags, 1);
% mixedJCVolume = zeros(nVolume, 1);
% parfor j = 1:nVolume
%     volumeData = load(volumeFile, volumeFields{j});
%     volumeDataFields = fieldnames(volumeData);
%     volumeMap = volumeData.(volumeDataFields{1});
%     mixedJCVolume(j, 1) = volumeMap.JC;
% end
% [JCVolume, sortJCIdx] = sort(mixedJCVolume, 'descend');
% clear mixedJCVolume
% sortedFields = volumeFields(sortJCIdx);
% countsVolume = zeros(nVolume, nSample^2);
% flagsVolume = zeros(nVolume, nSample^2);
% qVolume = zeros(6, nSample^2, nVolume);
% parfor j = 1:nVolume
%     volumeData = load(volumeFile, sortedFields{j});
%     volumeDataFields = fieldnames(volumeData);
%     volumeMap = volumeData.(volumeDataFields{1});
%     countsVolume(j,:) = reshape(volumeMap.counts, 1, nSample^2);
%     flagsVolume(j,:) = reshape(volumeMap.flags, 1, nSample^2);
%     qVolume(:,:,j) = volumeMap.q;
% end

%% Save JC Volume Data
% save(volumeDataFile, "JCVolume", "countsVolume", "flagsVolume", "qVolume", "-v7.3");

%% Load JC Volume Data
% load(volumeDataFile, "JCVolume", "flagsVolume", "qVolume");
% nVolume = length(JCVolume);
% nSample = sqrt(size(flagsVolume, 2));
% disp("Successfully loaded volume data from MAT file!")

%% JC Animation
% fig3 = figure('Position', [200 100 1200 750]);
% 
% ax3 = gca;
% axis(ax3, 'equal');
% grid(ax3, 'off');
% hold(ax3, 'on');
% 
% xlabel(ax3, "$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel(ax3, "$y$ [E-M ndim]", 'Interpreter', 'latex')
% view(ax3, 2)
% 
% % [Earth, ~, ~, ~] = plot3DBodyVid(ax3, "Earth", RE/lstarEM, [-muEM, 0, 0]);
% % set(Earth, 'DisplayName', "Earth")
% [Moon, ~, ~, ~] = plot3DBodyVid(ax3, "Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% 
% % scatSize = 1.75;
% scatSize = 7.5;
% hist3 = scatter(NaN(1, nSample^2), NaN(1, nSample^2), scatSize, repmat(colorMap(7,:), nSample^2, 1), 'filled', 'HandleVisibility', 'off');
% 
% % xlim(ax3, [-1.25 1.25])
% % ylim(ax3, [-1.25 1.25])
% xlim(ax3, [1-muEM-0.3 1-muEM+0.3])
% ylim(ax3, [-0.3 0.3])
% 
% scatter(nan, nan, 20, colorMap(7,:), 'filled', 'DisplayName', "Capture")
% scatter(nan, nan, 20, colorMap(8,:), 'filled', 'DisplayName', "Impact")
% scatter(nan, nan, 20, colorMap(10,:), 'filled', 'DisplayName', "ZVC")
% colormap(colorMap(1:6,:))
% cb3 = colorbar;
% clim([-0.5 5.5])
% cb3.Ticks = 0:5;
% cb3.TickLabels = [string(0:4), "5+"];
% ylabel(cb3, "Periapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb3.Label.Position = cb3.Label.Position+[-2.3 3.1 0];
% leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg3.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% 
% tit3 = annotation(fig3, 'textbox', [0.4 0.93 0.2 0.06], 'String', sprintf('JC = %.4f', JCVolume(1)), 'FontName', 'Times New Roman', 'FontSize', 18, 'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% 
% % v = VideoWriter('EscapeCR3BP_JCVolume_Moon_peri_retro.mp4', 'MPEG-4');
% % v.FrameRate = 40;
% % v.Quality = 100;
% % open(v);
% 
% for j = 1:nVolume
%     set(hist3, 'XData', qVolume(1,:,j), 'YData', qVolume(2,:,j), 'CData', colorMap(flagsVolume(j,:)+1,:));
% 
%     set(tit3, 'String', sprintf('JC = %.4f', JCVolume(j)));
% 
%     drawnow limitrate;
%     pause(0.1)
%     % writeVideo(v, getframe(fig3));
% end
% 
% % writeVideo(v, getframe(fig3));
% % close(v);

%% Area Analysis
% esc0 = zeros(nVolume, 1);
% esc1 = zeros(nVolume, 1);
% esc2 = zeros(nVolume, 1);
% esc3 = zeros(nVolume, 1);
% esc4 = zeros(nVolume, 1);
% esc5 = zeros(nVolume, 1);
% capture = zeros(nVolume, 1);
% impact = zeros(nVolume, 1);
% parfor j = 1:nVolume
%     esc0(j,1) = length(find(flagsVolume(j,:) == 0));
%     esc1(j,1) = length(find(flagsVolume(j,:) == 1));
%     esc2(j,1) = length(find(flagsVolume(j,:) == 2));
%     esc3(j,1) = length(find(flagsVolume(j,:) == 3));
%     esc4(j,1) = length(find(flagsVolume(j,:) == 4));
%     esc5(j,1) = length(find(flagsVolume(j,:) == 5));
%     capture(j,1) = length(find(flagsVolume(j,:) == 6));
%     impact(j,1) = length(find(flagsVolume(j,:) == 7));
% end
% escape = esc0+esc1+esc2+esc3+esc4+esc5;

% fig4 = figure("Position", [200 100 1200 750]);
% hold on
% plot(JCVolume, 100.*esc0./nSample^2, 'Color', colorMap(1,:), 'DisplayName', "Direct Esc.")
% plot(JCVolume, 100.*esc1./nSample^2, 'Color', colorMap(2,:), 'DisplayName', "Esc. +1")
% plot(JCVolume, 100.*esc2./nSample^2, 'Color', colorMap(3,:), 'DisplayName', "Esc. +2")
% plot(JCVolume, 100.*esc3./nSample^2, 'Color', colorMap(4,:), 'DisplayName', "Esc. +3")
% plot(JCVolume, 100.*esc4./nSample^2, 'Color', colorMap(5,:), 'DisplayName', "Esc. +4")
% plot(JCVolume, 100.*esc5./nSample^2, 'Color', colorMap(6,:), 'DisplayName', "Esc. +$\geq$5")
% plot(JCVolume, 100.*escape./nSample^2, 'w:', 'DisplayName', "Total Esc.")
% % xlim([2.95 3.15])
% xlabel("JC [E-M ndim]", 'Interpreter', 'latex')
% ylabel("Area [\%]", 'Interpreter', 'latex')
% title("Area Comparison", 'Interpreter', 'latex')
% leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax4 = gca;
% ax4.SortMethod = 'childorder';
% % exportgraphics(fig4, 'EscapeCR3BP_4.png', 'BackgroundColor', 'k')

%% Import Assisted Escape Analysis Data
% assistedData = load('../PhDScripts/Output/AssistedEscapeAnalysisCR3BP.mat');
% Deltav0s = assistedData.Deltav0s;

%% Assisted Escape Figure
% DeltavColors = hot(1000);
% maxDeltav = 250;
% pointColors = zeros(length(Deltav0s), 3);
% parfor j = 1:n^2
%     if ~isnan(Deltav0s(j))
%         pointColors(j,:) = getColor(DeltavColors, Deltav0s(j)*1000*lstarEM/tstarEM, [0, maxDeltav]);
%     else
%         pointColors(j,:) = [1, 1, 1];
%     end
% end
% 
% fig6 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xGrid, yGrid, 1.75, pointColors, 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter(nan, nan, 1.75, 'w', 'filled', 'DisplayName', "Infeasible")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
% colormap(hot)
% cb6 = colorbar;
% clim([0 maxDeltav])
% cb6.Ticks = 0:50:maxDeltav;
% cb6.TickLabels = [string(0:50:maxDeltav-50), string(maxDeltav)+"+"];
% ylabel(cb6, "$\Delta v$ [m/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb6.Label.Position = cb6.Label.Position+[-3.5 130 0];
% leg6 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% leg6.Position = leg6.Position+[0.1 0 0 0];
% drawnow;
% set(leg6.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax6 = gca;
% ax6.SortMethod = 'childorder';
% % exportgraphics(fig6, 'EscapeCR3BP_6.png','BackgroundColor', 'k')

%% Test Trajectory
% xSample = 0.864228;
% ySample = -0.0576152;
% idx = find((abs(xGrid-xSample) < 1E-5) & (abs(yGrid-ySample) < 1E-5));
% qOrig = map.q(:,idx);
% disp("Original IC: ["+qOrig(1)+", "+qOrig(2)+", "+qOrig(3)+", "+qOrig(4)+", "+qOrig(5)+", "+qOrig(6)+"]")
% tauOrig = 12*pi;
% solOrig = ode89(odeCR3BPEM, [0 tauOrig], qOrig, odeOpts);
% vOrig = norm(qOrig(4:5));
% vhat = qOrig(4:5)./vOrig;
% qAssist = qOrig;
% qAssist(4:5) = (vOrig+Deltav0s(idx)).*vhat;
% disp("Delta-v: "+Deltav0s(idx)*1000*lstarEM/tstarEM+" m/s")
% disp("Assisted IC: ["+qAssist(1)+", "+qAssist(2)+", "+qAssist(3)+", "+qAssist(4)+", "+qAssist(5)+", "+qAssist(6)+"]")
% tauAssist = 6*pi;
% solAssist = ode89(odeCR3BPEM, [0 tauAssist], qAssist, odeOpts);
% 
% fig7 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% scatter3(solOrig.y(1,1), solOrig.y(2,1), solOrig.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
% p71 = plot3WithArrows(solOrig.y(1,:), solOrig.y(2,:), solOrig.y(3,:), 'Color', colorMap(flags(idx)+1,:));
% set(p71, 'DisplayName', "Original Traj.")
% p72 = plot3WithArrows(solAssist.y(1,:), solAssist.y(2,:), solAssist.y(3,:), 'g');
% set(p72, 'DisplayName', "Assisted Traj.")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% % axis([-3 3 -3 3])
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg7 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg7.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax7 = gca;
% ax7.SortMethod = 'childorder';
% % exportgraphics(fig7, 'EscapeCR3BP_7.png','BackgroundColor', 'k')