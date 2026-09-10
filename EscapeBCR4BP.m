%%% EscapeBCR4BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 19 August 2026
%%% U: 24 August 2026

clear

%% Import Map Data
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
% colorMap(10,:) = [0, 0, 0]; % ZVC
% colorMap(9,:) = [1, 1, 1]; % Invalid apse

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
% % title("Earth-Moon Rot.: JC = "+JC+" | $\theta_{S}$ = "+thetaS+"$^{\circ}$", 'Interpreter', 'latex')
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

%% Import Sun Angle Volume Data
% volumeFile = 'E:/ApseMapData/BCR4BPthetaVolume_1_peri_pro_500_3.0663.mat';
% volumeDataFile = 'BCR4BPthetaVolume_1_peri_pro_500_3.0663.mat';

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
% load(volumeDataFile, "thetaVolume", "flagsVolume", "qVolume");
% nVolume = length(thetaVolume);
% nSample = size(qVolume, 2);
% disp("Successfully loaded volume data from MAT file!")

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