%%% EscapeBCR4BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 19 August 2026
%%% U: 15 October 2026

clear

%% Import Map Data
mapsDataCRBP = load('../PhDScripts/Output/ApseMaps/CR3BP_1_peri_pro_500_3.0663.mat');
mapFieldsCR3BP = fieldnames(mapsDataCRBP);
mapCR3BP = mapsDataCRBP.(mapFieldsCR3BP{1});
flagsCR3BP = mapCR3BP.flags;

mapsData = load('../PhDScripts/Output/ApseMaps/BCR4BP_1_peri_pro_500_3.0663_0.0.mat');
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
thetaS = rad2deg(map.thetaS);
disp(primary+"-centered "+grade+" "+apse+" map: JC = "+JC+" | Sun angle = "+thetaS)

n = length(map.flags);
xGrid = map.q(1,:);
yGrid = map.q(2,:);
flags = map.flags;
periapses = map.periapses;
periapsesIndices = map.periapsesIndices;

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

m4 = 3.2890055981024733E5;
a4 = 3.888204672248271E2;

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
% colorMap(10,:) = [0, 0, 0]; % ZVC
% colorMap(9,:) = [1, 1, 1]; % Invalid apse

%% CR3BP Map
% fig0 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xGrid, yGrid, 1.75, colorMap(flagsCR3BP+1,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter(nan, nan, 20, colorMap(7,:), 'filled', 'DisplayName', "Capture")
% scatter(nan, nan, 20, colorMap(8,:), 'filled', 'DisplayName', "Impact")
% scatter(nan, nan, 20, colorMap(10,:), 'filled', 'DisplayName', "ZVC")
% axis equal
% if primary == "Moon"
%     axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% else
%     axis([-1.25 1.25 -1.25 1.25])
% end
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
% colormap(colorMap(1:6,:))
% cb0 = colorbar;
% clim([-0.5 5.5])
% cb0.Ticks = 0:5;
% cb0.TickLabels = [string(0:4), "5+"];
% ylabel(cb0, "Periapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% % ylabel(cb0, "Apoapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb0.Label.Position = cb0.Label.Position+[-2.3 3.1 0];
% leg0 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg0.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax0 = gca;
% ax0.SortMethod = 'childorder';

%% Map
% fig1 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xGrid, yGrid, 1.75, colorMap(flags+1,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter(nan, nan, 20, colorMap(7,:), 'filled', 'DisplayName', "Capture")
% scatter(nan, nan, 20, colorMap(8,:), 'filled', 'DisplayName', "Impact")
% scatter(nan, nan, 20, colorMap(10,:), 'filled', 'DisplayName', "ZVC")
% axis equal
% if primary == "Moon"
%     axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% else
%     axis([-1.25 1.25 -1.25 1.25])
% end
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC+" | $\theta_{S}$ = "+thetaS+"$^{\circ}$", 'Interpreter', 'latex')
% colormap(colorMap(1:6,:))
% cb1 = colorbar;
% clim([-0.5 5.5])
% cb1.Ticks = 0:5;
% cb1.TickLabels = [string(0:4), "5+"];
% ylabel(cb1, "Periapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% % ylabel(cb1, "Apoapses", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb1.Label.Position = cb1.Label.Position+[-2.3 3.1 0];
% leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg1.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax1 = gca;
% ax1.SortMethod = 'childorder';
% % exportgraphics(fig1, 'EscapeBCR4BP_1.png', 'BackgroundColor', 'k')
% % exportgraphics(fig1, 'EscapeBCR4BP_1.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Difference Colormap
difColorMap(1,:) = [0.267, 0.005, 0.329]; % Robust escape
difColorMap(2,:) = [0, 0.5, 0.5]; % Solar-inhibited escape
difColorMap(3,:) = [0.5, 0.5, 0]; % Solar-assisted escape
difColorMap(4,:) = [0.78, 0.72, 0.66]; % Robust capture
difColorMap(5,:) = [0, 0, 0]; % Invalid apse
difColorMap(6,:) = [1, 1, 1]; % ZVC

%% Difference Map
% flagsDif = flags;
% flagsDif((flags < 6) & (flagsCR3BP < 6)) = 1; % Robust escape
% flagsDif((flags >= 6) & (flagsCR3BP < 6)) = 2; % Solar-inhibited escape
% flagsDif((flags < 6) & (flagsCR3BP >= 6)) = 3; % Solar-assisted escape
% flagsDif((flags >= 6) & (flagsCR3BP >= 6)) = 4; % Robust capture
% flagsDif(flags == 8) = 5; % Invalid apse
% flagsDif(flags == 9) = 6; % ZVC
% 
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xGrid, yGrid, 1.75, difColorMap(flagsDif,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter(nan, nan, 20, difColorMap(1,:), 'filled', 'DisplayName', "Robust Escape")
% scatter(nan, nan, 20, difColorMap(2,:), 'filled', 'DisplayName', "Solar-Inhibited Escape")
% scatter(nan, nan, 20, difColorMap(3,:), 'filled', 'DisplayName', "Solar-Assisted Escape")
% scatter(nan, nan, 20, difColorMap(4,:), 'filled', 'DisplayName', "Robust Capture/Impact")
% scatter(nan, nan, 20, difColorMap(6,:), 'filled', 'DisplayName', "ZVC")
% axis equal
% if primary == "Moon"
%     axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% else
%     axis([-1.25 1.25 -1.25 1.25])
% end
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC+" | $\theta_{S}$ = "+thetaS+"$^{\circ}$", 'Interpreter', 'latex')
% leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg2.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax2 = gca;
% ax2.SortMethod = 'childorder';
% % exportgraphics(fig2, 'EscapeBCR4BP_2.png', 'BackgroundColor', 'k')
% % exportgraphics(fig2, 'EscapeBCR4BP_2.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Propagators
odeCR3BPEM = @(t,r) ODE_CR3BP_JC(t, r, muEM, m4, a4, mS, mstarEM);
odeBCR4BPEM = @(t,r) ODE_BCR4BP(t, r, muEM, m4, a4, mS, mstarEM);
odeOpts = odeset('RelTol', 1E-12, 'AbsTol', 1E-12);

%% Test Trajectory
% xSample = -0.0216433;
% ySample = 0.502605;
% idx = find((abs(xGrid-xSample) < 1E-5) & (abs(yGrid-ySample) < 1E-5))
% % idx = 120003; % Solar-Assisted Escape
% 
% q = [map.q(:,idx); JC];
% disp("Sample IC: ["+q(1)+", "+q(2)+", "+q(3)+", "+q(4)+", "+q(5)+", "+q(6)+", "+q(7)+"]")
% tau = 12*pi;
% solCR3BP = ode89(odeCR3BPEM, [0 tau], q, odeOpts);
% solBCR4BP = ode89(odeBCR4BPEM, [0, tau], q(1:7), odeOpts);
% instJC = zeros(size(solBCR4BP.x));
% for j = 1:length(instJC)
%     instJC(j) = getJC(muEM, solBCR4BP.y(:,j));
% end
% periIndices = find(abs(periapsesIndices-idx) < 1E-5);
% peris = periapses(:,periIndices(1):periIndices(end));
% thetaSdot = sqrt((mS/mstarEM+1)/(a4^3))-1;
% periTimes = peris(7,:)./thetaSdot;
% periJCs = zeros(size(periTimes));
% for j = 1:length(periTimes)
%     periJCs(j) = getJC(muEM, peris(:,j));
% end

% fig4 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% scatter3(solCR3BP.y(1,1), solCR3BP.y(2,1), solCR3BP.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
% p41 = plot3WithArrows(solCR3BP.y(1,:), solCR3BP.y(2,:), solCR3BP.y(3,:), 'Color', colorMap(flagsCR3BP(idx)+1,:), 'NumArrows', 5);
% set(p41, 'DisplayName', "CR3BP")
% p42 = plot3WithArrows(solBCR4BP.y(1,:), solBCR4BP.y(2,:), solBCR4BP.y(3,:), 'Color', colorMap(flags(idx)+1,:), 'NumArrows', 10, 'ArrowScale', 0.2);
% set(p42, 'DisplayName', "$\theta_{S}=0$")
% scatter3(peris(1,:), peris(2,:), peris(3,:), 50, 'm', 'filled', 'DisplayName', "Periapses")
% axis equal
% if primary == "Moon"
%     axis([1-muEM-0.2 1-muEM+0.2 -0.2 0.2])
% else
%     axis([-1.25 1.25 -1.25 1.25])
% end
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg4.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax4 = gca;
% ax4.SortMethod = 'childorder';
% % exportgraphics(fig4, 'EscapeBCR4BP_4.png','BackgroundColor', 'k')

% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% yline(JC, 'w--', 'LineWidth', 1, 'DisplayName', "CR3BP");
% plot(solBCR4BP.x, instJC, 'DisplayName', "Inst. JC")
% plot(solCR3BP.x, solCR3BP.y(8,:), 'y', 'DisplayName', "Est. JC")
% xline(2*pi/-thetaSdot, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
% xline(4*pi/-thetaSdot, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
% xline(6*pi/-thetaSdot, 'w:', 'LineWidth', 1, 'HandleVisibility', 'off')
% scatter(periTimes, periJCs, 50, 'm', 'filled', 'DisplayName', "Periapses")
% axis([0 25 JC-0.025 JC+0.025])
% xlabel("$\tau$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("JC", 'Interpreter', 'latex')
% leg5 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax5 = gca;
% ax5.SortMethod = 'childorder';
% % exportgraphics(fig5, 'EscapeBCR4BP_5.png', 'BackgroundColor', 'k')

%% Import Sun Angle Volume Data
volumeFile = 'E:/ApseMapData/BCR4BPthetaVolume_1_peri_pro_500_3.0663.mat';
volumeDataFile = 'BCR4BPthetaVolume_1_peri_pro_500_3.0663_new.mat';

% volumeFields = who('-file', volumeFile);
% nVolume = length(volumeFields);
% sampleData = load(volumeFile, volumeFields{1});
% sampleFields = fieldnames(sampleData);
% sampleMap = sampleData.(sampleFields{1});
% rSample = size(sampleMap.flags, 1);
% cSample = size(sampleMap.flags, 2);
% mixedthetaVolume = zeros(nVolume, 1);
% parfor j = 1:nVolume
%     volumeData = load(volumeFile, volumeFields{j});
%     volumeDataFields = fieldnames(volumeData);
%     volumeMap = volumeData.(volumeDataFields{1});
%     mixedthetaVolume(j, 1) = volumeMap.thetaS;
% end
% [thetaVolume, sortthetaIdx] = sort(mixedthetaVolume, 'ascend');
% clear mixedthetaVolume
% sortedFields = volumeFields(sortthetaIdx);
% countsVolume = zeros(nVolume, rSample, cSample);
% flagsVolume = zeros(nVolume, rSample, cSample);
% qVolume = zeros(7, rSample*cSample, nVolume);
% parfor j = 1:nVolume
%     volumeData = load(volumeFile, sortedFields{j});
%     volumeDataFields = fieldnames(volumeData);
%     volumeMap = volumeData.(volumeDataFields{1});
%     countsVolume(j,:) = reshape(volumeMap.counts, 1, rSample, cSample);
%     flagsVolume(j,:) = reshape(volumeMap.flags, 1, rSample, cSample);
%     qVolume(:,:,j) = volumeMap.q;
% end

%% Save Sun Angle Volume Data
% save(volumeDataFile, "thetaVolume", "countsVolume", "flagsVolume", "qVolume", "-v7.3");

%% Load Sun Angle Volume Data
load(volumeDataFile, "thetaVolume", "flagsVolume", "qVolume");
nVolume = length(thetaVolume);
nSample = size(qVolume, 2);
disp("Successfully loaded volume data from MAT file!")

%% Sun Angle Animation
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
% [Earth, ~, ~, ~] = plot3DBodyVid(ax3, "Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% [Moon, ~, ~, ~] = plot3DBodyVid(ax3, "Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% 
% scatSize = 1.75;
% % scatSize = 7.5;
% hist3 = scatter(NaN(1, nSample), NaN(1, nSample), scatSize, repmat(colorMap(7,:), nSample, 1), 'filled', 'HandleVisibility', 'off');
% 
% xlim(ax3, [-1.25 1.25])
% ylim(ax3, [-1.25 1.25])
% % xlim(ax3, [1-muEM-0.3 1-muEM+0.3])
% % ylim(ax3, [-0.3 0.3])
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
% drawnow;
% set(leg3.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% 
% tit3 = annotation(fig3, 'textbox', [0.4 0.93 0.2 0.06], 'String', "$\theta_{S}="+rad2deg(thetaVolume(1))+"^{\circ}$", 'FontName', 'Times New Roman', 'FontSize', 18, 'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex');
% 
% % v = VideoWriter('EscapeBCR4BP_thetaVolume_Earth_peri_pro.mp4', 'MPEG-4');
% % v.FrameRate = 40;
% % v.Quality = 100;
% % open(v);
% 
% for j = 1:nVolume
%     set(hist3, 'XData', qVolume(1,:,j), 'YData', qVolume(2,:,j), 'CData', colorMap(flagsVolume(j,:)+1,:));
% 
%     set(tit3, 'String', "$\theta_{S}="+rad2deg(thetaVolume(j))+"^{\circ}$");
% 
%     drawnow limitrate;
%     pause(0.1)
%     % writeVideo(v, getframe(fig3));
% end
% 
% % writeVideo(v, getframe(fig3));
% % close(v);

%% Full Range Difference Map
flagsDifFull = flagsCR3BP;
flagsDifFull(all(flagsVolume < 6, 1).' & (flagsCR3BP < 6)) = 1; % Robust escape
flagsDifFull(any(flagsVolume >= 6, 1).' & (flagsCR3BP < 6)) = 2; % Solar-inhibited escape
flagsDifFull(any(flagsVolume < 6, 1).' & (flagsCR3BP >= 6)) = 3; % Solar-assisted escape
flagsDifFull(all(flagsVolume >= 6, 1).' & (flagsCR3BP >= 6)) = 4; % Robust capture
flagsDifFull(flagsCR3BP == 8) = 5; % Invalid apse
flagsDifFull(flagsCR3BP == 9) = 6; % ZVC

% fig6 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xGrid, yGrid, 1.75, difColorMap(flagsDifFull,:), 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter(nan, nan, 20, difColorMap(1,:), 'filled', 'DisplayName', "Robust Escape")
% scatter(nan, nan, 20, difColorMap(2,:), 'filled', 'DisplayName', "Solar-Inhibited Escape")
% scatter(nan, nan, 20, difColorMap(3,:), 'filled', 'DisplayName', "Solar-Assisted Escape")
% scatter(nan, nan, 20, difColorMap(4,:), 'filled', 'DisplayName', "Robust Capture/Impact")
% scatter(nan, nan, 20, difColorMap(6,:), 'filled', 'DisplayName', "ZVC")
% axis equal
% if primary == "Moon"
%     axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% else
%     axis([-1.25 1.25 -1.25 1.25])
% end
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
% leg6 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg6.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax6 = gca;
% ax6.SortMethod = 'childorder';
% % exportgraphics(fig6, 'EscapeBCR4BP_6.png', 'BackgroundColor', 'k')
% % exportgraphics(fig6, 'EscapeBCR4BP_6.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Area Analysis
esc0 = zeros(nVolume, 1);
esc1 = zeros(nVolume, 1);
esc2 = zeros(nVolume, 1);
esc3 = zeros(nVolume, 1);
esc4 = zeros(nVolume, 1);
esc5 = zeros(nVolume, 1);
capture = zeros(nVolume, 1);
impact = zeros(nVolume, 1);
parfor j = 1:nVolume
    esc0(j,1) = length(find(flagsVolume(j,:) == 0));
    esc1(j,1) = length(find(flagsVolume(j,:) == 1));
    esc2(j,1) = length(find(flagsVolume(j,:) == 2));
    esc3(j,1) = length(find(flagsVolume(j,:) == 3));
    esc4(j,1) = length(find(flagsVolume(j,:) == 4));
    esc5(j,1) = length(find(flagsVolume(j,:) == 5));
    capture(j,1) = length(find(flagsVolume(j,:) == 6));
    impact(j,1) = length(find(flagsVolume(j,:) == 7));
end
escape = esc0+esc1+esc2+esc3+esc4+esc5;
escCR3BP = length(find(flagsCR3BP < 6));
escRobust = length(find(flagsDifFull == 1));
escBCR4BP = length(find(flagsDifFull < 4));

fig7 = figure("Position", [200 100 1200 750]);
hold on
plot(thetaVolume.*180./pi, 100.*escape./nSample, 'b', 'DisplayName', "Total Esc.")
plot(thetaVolume.*180./pi, 100.*esc0./nSample, 'Color', colorMap(1,:), 'DisplayName', "Direct Esc.")
yline(100*escCR3BP/nSample, 'y--', 'DisplayName', "CR3BP Esc.")
yline(100*escRobust/nSample, 'w--', 'DisplayName', "Robust Esc.")
yline(100*escBCR4BP/nSample, 'g--', 'DisplayName', "BCR4BP Feas. Esc.")
xlim([0 360])
xlabel("$\theta_{S}$ [E-M ndim]", 'Interpreter', 'latex')
ylabel("Area [\%]", 'Interpreter', 'latex')
title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
leg7 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
ax7 = gca;
ax7.SortMethod = 'childorder';
% exportgraphics(fig7, 'EscapeBCR4BP_7.png', 'BackgroundColor', 'k')