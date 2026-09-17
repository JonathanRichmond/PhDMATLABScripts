%%% EscapeCR3BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 16 June 2026
%%% U: 17 September 2026

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
rHill = 3.8897077389538994; % Earth Hill sphere radius

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

%% Import Apse Data
% apseData = load('../PhDScripts/Output/ApseStatesCR3BP.mat');
% qOrbit = apseData.qOrbit;
% tOrbit = apseData.orbitP;
% initialStates = apseData.initialStates;
% periStates = apseData.periStates;
% periTimes = apseData.periTimes;
% xPeri = periStates(1,:)';
% yPeri = periStates(2,:)';
% 
% trajEscE = apseData.escE;
% trajEscv = apseData.escv;
% trajq0 = apseData.q0;
% trajt1 = apseData.t1;
% trajDeltav1 = apseData.Deltav1;
% trajq1 = apseData.q1;
% trajt2 = apseData.t2;
% trajDeltav2 = apseData.Deltav2;
% trajq2 = apseData.q2;
% trajNewJC = apseData.newJC;
% trajNewEscE = apseData.newEscE;
% trajNewEscv = apseData.newEscv;
% trajPeriStates = apseData.peris;
% trajApoStates = apseData.apos;

% % xTrajSample = -0.0650211;
% % yTrajSample = 0.208377;
% % trajIdx = find((abs(xPeri-xTrajSample) < 1E-5) & (abs(yPeri-yTrajSample) < 1E-5))
% % trajIdx = 63 % Indirect
% trajIdx = 99 % Capture
% % trajIdx = 22 % Capture with option for high-energy maneuver

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
% % scatter(xPeri, yPeri, 30, 'm', 'filled', 'DisplayName', "Manifold Peris.")
% % scatter(xPeri(trajIdx), yPeri(trajIdx), 50, 'g', 'filled', 'DisplayName', "Sample")
% % scatter(sol1.y(1,1), sol1.y(2,1), 50, 'c', 'filled', 'HandleVisibility', 'off')
% % scatter(sol2.y(1,1), sol2.y(2,1), 50, 'c', 'filled', 'HandleVisibility', 'off')
% % scatter3(trajPeriStates(1,:), trajPeriStates(2,:), trajPeriStates(3,:), 50, 'c', 'o', 'filled', 'HandleVisibility', 'off')
% axis equal
% if primary == "Moon"
%     axis([1-muEM-0.2 1-muEM+0.2 -0.2 0.2])
% else
%     axis([-1.25 1.25 -1.25 1.25])
% end
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
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
% % exportgraphics(fig1, 'EscapeCR3BP_1.png', 'BackgroundColor', 'k')
% % exportgraphics(fig1, 'EscapeCR3BP_1.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Propagators
odeCR3BPEM = @(t,r) ODE_CR3BP(t, r, muEM);
odeOpts = odeset('RelTol', 1E-12, 'AbsTol', 1E-12);

%% Test Trajectory
% xSample = -0.158317;
% ySample = 0.102204;
% idx = find((abs(xGrid-xSample) < 1E-5) & (abs(yGrid-ySample) < 1E-5))
% idx = 46420 % Direct
% idx = 124010 % Indirect
% idx = 133530 % Failure

% q = map.q(:,idx);
% % q = initialStates(:,trajIdx);
% disp("Sample IC: ["+q(1)+", "+q(2)+", "+q(3)+", "+q(4)+", "+q(5)+", "+q(6)+"]")
% tau = 2*pi;
% sol = ode89(odeCR3BPEM, [0 tau], q, odeOpts);

% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% scatter3(sol.y(1,1), sol.y(2,1), sol.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
% p21 = plot3WithArrows(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', colorMap(flags(idx)+1,:));
% % p21 = plot3WithArrows(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', 'm');
% if flags(idx) > 0
%     periIndices = find(abs(periapsesIndices-idx) < 1E-5);
%     peris = periapses(:,periIndices(1):periIndices(end));
%     scatter3(peris(1,:), peris(2,:), peris(3,:), 50, 'm', 'filled', 'DisplayName', "Periapses")
% end
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
% set(gca, 'Color', 'w');
% view(2)
% hold off
% ax2 = gca;
% ax2.SortMethod = 'childorder';
% % exportgraphics(fig2, 'EscapeCR3BP_2.png','BackgroundColor', 'k')

%% Trajectory Analysis
% solOrbit = ode89(odeCR3BPEM, [0 tOrbit], qOrbit, odeOpts);
% disp("Initial State: ["+trajq0(1)+", "+trajq0(2)+", "+trajq0(3)+", "+trajq0(4)+", "+trajq0(5)+", "+trajq0(6)+"]")
% disp("Original departure energy:"+trajEscE)
% disp("Original departure velocity:"+trajEscv)
% tau = 12*pi;
% sol = ode89(odeCR3BPEM, [0 periTimes(trajIdx)+tau], initialStates(:,trajIdx), odeOpts);
% sol0 = ode89(odeCR3BPEM, [0 periTimes(trajIdx)+trajt1], initialStates(:,trajIdx), odeOpts);
% disp("Delta-v 1: "+trajDeltav1*1000*lstarEM/tstarEM+" m/s")
% disp("Indirect State: ["+trajq1(1)+", "+trajq1(2)+", "+trajq1(3)+", "+trajq1(4)+", "+trajq1(5)+", "+trajq1(6)+"]")
% if trajt2 > 0
%     sol1 = ode89(odeCR3BPEM, [trajt1 trajt1+trajt2], trajq1, odeOpts);
% end
% sol15 = ode89(odeCR3BPEM, [trajt1 trajt1+trajt2+1.5*pi], trajq1, odeOpts);
% disp("Delta-v 2: "+trajDeltav2*1000*lstarEM/tstarEM+" m/s")
% disp("Escape State: ["+trajq2(1)+", "+trajq2(2)+", "+trajq2(3)+", "+trajq2(4)+", "+trajq2(5)+", "+trajq2(6)+"]")
% sol2 = ode89(odeCR3BPEM, [trajt1+trajt2 trajt1+trajt2+3*pi], trajq2, odeOpts);
% disp("New JC: "+trajNewJC)
% disp("New departure energy:"+trajNewEscE)
% disp("New departure velocity:"+trajNewEscv)
% disp("Departure metric:"+(trajNewEscv-trajEscv)/((abs(trajDeltav1)+abs(trajDeltav2))*lstarEM/tstarEM))

% fig13 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% plot3(rHill.*cos(linspace(0, 2*pi, 1001)), rHill.*sin(linspace(0, 2*pi, 1001)), zeros(1,1001), 'w:', 'DisplayName', "Hills Sphere")
% % plot3(solOrbit.y(1,:), solOrbit.y(2,:), solOrbit.y(3,:), 'g:', 'DisplayName', "Orbit")
% p131 = plot3WithArrows(solOrbit.y(1,:), solOrbit.y(2,:), solOrbit.y(3,:), 'LineType', ':', 'Color', 'g', 'NumArrows', 2, 'ArrowScale', 5);
% set(p131, 'DisplayName', "Orbit")
% % plot3(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'm', 'DisplayName', "Orig. Traj.")
% p132 = plot3WithArrows(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', 'm', 'NumArrows', 10, 'ArrowScale', 1);
% set(p132, 'DisplayName', "Orig. Traj.")
% % scatter(sol1.y(1,1), sol1.y(2,1), 50, 'c', 'filled', 'HandleVisibility', 'off')
% % scatter3(trajPeriStates(1,:), trajPeriStates(2,:), trajPeriStates(3,:), 50, 'c', 'o', 'filled', 'HandleVisibility', 'off')
% % scatter3(trajq0(1), trajq0(2), trajq0(3), 50, 'g', 'o', 'filled', 'DisplayName', "Sample Perigee")
% plot3(sol0.y(1,:), sol0.y(2,:), sol0.y(3,:), 'b', 'DisplayName', "New Traj.")
% % % plot3(sol15.y(1,:), sol15.y(2,:), sol15.y(3,:), 'b', 'HandleVisibility', 'off')
% % p134 = plot3WithArrows(sol15.y(1,:), sol15.y(2,:), sol15.y(3,:), 'Color', 'b', 'NumArrows', 10, 'ArrowScale', 0.5);
% % set(p134, 'HandleVisibility', 'off')
% % scatter(sol2.y(1,1), sol2.y(2,1), 50, 'c', 'filled', 'HandleVisibility', 'off')
% if trajt2 > 0
%     % plot3(sol1.y(1,:), sol1.y(2,:), sol1.y(3,:), 'b', 'HandleVisibility', 'off')
%     p134 = plot3WithArrows(sol1.y(1,:), sol1.y(2,:), sol1.y(3,:), 'Color', 'b', 'NumArrows', 5, 'ArrowScale', 1.2);
%     set(p134, 'HandleVisibility', 'off')
% end
% % plot3(sol2.y(1,:), sol2.y(2,:), sol2.y(3,:), 'b', 'HandleVisibility', 'off')
% p135 = plot3WithArrows(sol2.y(1,:), sol2.y(2,:), sol2.y(3,:), 'Color', 'b', 'NumArrows', 2, 'ArrowScale', 1);
% set(p135, 'HandleVisibility', 'off')
% if abs(trajDeltav1) > 0
%     scatter3(sol1.y(1,1), sol1.y(2,1), sol1.y(3,1), 50, 'w', 's', 'filled', 'DisplayName', "Man. 1")
% end
% if abs(trajDeltav2) > 0
%     scatter3(sol2.y(1,1), sol2.y(2,1), sol2.y(3,1), 50, 'w', '^', 'filled', 'DisplayName', "Man. 2")
% end
% axis equal
% % axis([-1.25 1.25 -1.25 1.25])
% % axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg13 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg13.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax13 = gca;
% ax13.SortMethod = 'childorder';
% % exportgraphics(fig13, 'EscapeCR3BP_13.png','BackgroundColor', 'k')
% % exportgraphics(fig13, 'EscapeCR3BP_13.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% High-Energy Analysis
% trajApo = trajApoStates(:,9);
% disp("Initial State: ["+trajApo(1)+", "+trajApo(2)+", "+trajApo(3)+", "+trajApo(4)+", "+trajApo(5)+", "+trajApo(6)+"]")
% tau = 12*pi;
% sol = ode89(odeCR3BPEM, [0 tau], trajApo, odeOpts);
% trajPeri = trajPeriStates(:,10);
% vMag = norm(trajPeri(4:6));
% vhat = trajPeri(4:6)./vMag;
% DeltavHigh = 1000*tstarEM/lstarEM/1000;
% qHigh1 = trajPeri;
% qHigh1(4:6) = vhat.*(vMag+DeltavHigh);
% disp("Delta-v 1: "+DeltavHigh*1000*lstarEM/tstarEM+" m/s")
% disp("Altered Perigee State: ["+qHigh1(1)+", "+qHigh1(2)+", "+qHigh1(3)+", "+qHigh1(4)+", "+qHigh1(5)+", "+qHigh1(6)+"]")
% tau1 = 2*pi;
% sol1 = ode89(odeCR3BPEM, [0 tau1], qHigh1, odeOpts);
% vMagApo = norm(trajApo(4:6));
% vhatApo = trajApo(4:6)./vMagApo;
% DeltavApo = 500*tstarEM/lstarEM/1000;
% qApo = trajApo;
% qApo(4:6) = vhatApo.*(vMagApo+DeltavApo);
% disp("Delta-v Apo: "+DeltavApo*1000*lstarEM/tstarEM+" m/s")
% disp("Altered Apogee State: ["+qApo(1)+", "+qApo(2)+", "+qApo(3)+", "+qApo(4)+", "+qApo(5)+", "+qApo(6)+"]")
% tauApo = 0.5*pi;
% solApo = ode89(odeCR3BPEM, [0 tauApo], qApo, odeOpts);
% newIdx = 45;
% trajPeriNew = solApo.y(:,newIdx);
% vMagNew = norm(trajPeriNew(4:6));
% vhatNew = trajPeriNew(4:6)./vMagNew;
% DeltavHighNew = DeltavHigh-DeltavApo;
% qHigh2 = trajPeriNew;
% qHigh2(4:6) = vhatNew.*(vMagNew+DeltavHighNew);
% disp("Delta-v 2: "+DeltavHighNew*1000*lstarEM/tstarEM+" m/s")
% disp("New Altered Perigee State: ["+qHigh2(1)+", "+qHigh2(2)+", "+qHigh2(3)+", "+qHigh2(4)+", "+qHigh2(5)+", "+qHigh2(6)+"]")
% tau2 = 2*pi;
% sol2 = ode89(odeCR3BPEM, [0 tau2], qHigh2, odeOpts);
% 
% fig14 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% % plot3(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'm', 'DisplayName', "Orig. Traj.")
% p141 = plot3WithArrows(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', 'm', 'NumArrows', 10, 'ArrowScale', 1);
% set(p141, 'DisplayName', "Orig. Traj.")
% scatter3(trajApo(1), trajApo(2), trajApo(3), 50, 'r', 'o', 'filled', 'DisplayName', "IC")
% scatter3(trajPeri(1), trajPeri(2), trajPeri(3), 50, 'g', 'o', 'filled', 'DisplayName', "Sample Peri.")
% % plot3(sol1.y(1,:), sol1.y(2,:), sol1.y(3,:), 'b', 'DisplayName', "1-Man. Traj.")
% p142 = plot3WithArrows(sol1.y(1,:), sol1.y(2,:), sol1.y(3,:), 'Color', 'b', 'NumArrows', 3, 'ArrowScale', 0.1);
% set(p142, 'DisplayName', "1-Man. Traj.")
% scatter3(sol1.y(1,1), sol1.y(2,1), sol1.y(3,1), 50, 'w', '^', 'filled', 'DisplayName', "Peri. Man.")
% % plot3(solApo.y(1,:), solApo.y(2,:), solApo.y(3,:), 'g', 'DisplayName', "2-Man. Traj.")
% p143 = plot3WithArrows(solApo.y(1,1:newIdx), solApo.y(2,1:newIdx), solApo.y(3,1:newIdx), 'Color', 'g', 'NumArrows', 1, 'ArrowScale', 3);
% set(p143, 'DisplayName', "2-Man. Traj.")
% % plot3(sol2.y(1,:), sol2.y(2,:), sol2.y(3,:), 'g', 'HandleVisibility', 'off')
% p144 = plot3WithArrows(sol2.y(1,:), sol2.y(2,:), sol2.y(3,:), 'Color', 'g', 'NumArrows', 2, 'ArrowScale', 0.1);
% set(p144, 'HandleVisibility', 'off')
% scatter3(solApo.y(1,1), solApo.y(2,1), solApo.y(3,1), 50, 'w', 's', 'filled', 'DisplayName', "Apo. Man.")
% scatter3(sol2.y(1,1), sol2.y(2,1), sol2.y(3,1), 50, 'w', '*', 'DisplayName', "New Peri. Man.")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% % axis([1-muEM-0.3 1-muEM+0.3 -0.3 0.3])
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg14 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg14.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax14 = gca;
% ax14.SortMethod = 'childorder';
% % exportgraphics(fig14, 'EscapeCR3BP_14.png','BackgroundColor', 'k')

%% Import Escape Analysis Data
% analysisData = load('../PhDScripts/Output/EscapeAnalysisCR3BP.mat');
% q0s_filt0 = analysisData.esc0q0;
% tfs_filt0 = analysisData.esc0tf;
% Es_filt0 = analysisData.esc0E;
% ds_filt0 = analysisData.esc0d;
% nFilt0 = length(Es_filt0);
% q0s_filt1 = analysisData.esc1q0;
% tfs_filt1 = analysisData.esc1tf;
% Es_filt1 = analysisData.esc1E;
% nFilt1 = length(Es_filt1);
% E_min = min([Es_filt0; Es_filt1]);
% E_max = max([Es_filt0; Es_filt1]);
% d_min = min(ds_filt0);
% d_max = max(ds_filt0);
% 
% Deltav2s = analysisData.Deltav2s;
% Deltav2bs = analysisData.Deltav2bs;
% escvs = analysisData.Escapevs;
% escvbs = analysisData.Escapevbs;
% % DeltaEs = analysisData.DeltaEs;
% gammas = analysisData.metrics;
% gammabs = analysisData.metricbs;
% 
% escE = analysisData.escE;
% escv = analysisData.escv;
% Deltav1 = analysisData.Deltav1;
% Deltav2 = analysisData.Deltav2;
% qMan = analysisData.qMan;
% newJC = analysisData.newJC;
% newEscE = analysisData.newEscE;
% newEscv = analysisData.newEscv;

%% Escape Analysis
% fig8 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(Deltav2s.*1000.*lstarEM./tstarEM, escvs, 20, 'filled', 'DisplayName', "Direct")
% % scatter(Deltav2bs.*1000.*lstarEM./tstarEM, escvbs, 20, 'filled', 'DisplayName', "Indirect")
% % scatter(Deltav2s.*1000.*lstarEM./tstarEM, gammas, 20, 'filled', 'DisplayName', "Direct")
% % scatter(Deltav2bs.*1000.*lstarEM./tstarEM, gammabs, 20, 'filled', 'DisplayName', "Indirect")
% % scatter(Deltav2s.*1000.*lstarEM./tstarEM, DeltaEs, 20, 'filled', 'HandleVisibility', 'off')
% % scatter(Deltav2s.*1000.*lstarEM./tstarEM, flybys, 20, 'filled', 'HandleVisibility', 'off')
% % ylim([-60 60])
% xlabel("$\Delta v_{2}$ [m/s]", 'Interpreter', 'latex')
% % xlabel("JC", 'Interpreter', 'latex')
% % ylabel("$\mathcal{v}_{dep}$ [km/s]", 'Interpreter', 'latex')
% ylabel("$\gamma_{esc}$", 'Interpreter', 'latex')
% % ylabel("$d\mathcal{E}_{esc}/d\alpha$ [km/s]", 'Interpreter', 'latex')
% % ylabel("$r_{p}$ [km]", 'Interpreter', 'latex')
% title("Maneuver Optimization", 'Interpreter', 'latex')
% leg8 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax8 = gca;
% ax8.SortMethod = 'childorder';
% % exportgraphics(fig8, 'EscapeCR3BP_8.png', 'BackgroundColor', 'k')

% qOrig = map.q(:,idx);
% disp("IC: ["+qOrig(1)+", "+qOrig(2)+", "+qOrig(3)+", "+qOrig(4)+", "+qOrig(5)+", "+qOrig(6)+"]")
% disp("Original departure energy:"+escE)
% disp("Original departure velocity:"+escv)
% tauOrig = 12*pi;
% solOrig = ode89(odeCR3BPEM, [0 tauOrig], qOrig, odeOpts);

% % % Deltav1
% % q1 = qMan;
% % disp("Delta-v 1: "+Deltav1*1000*lstarEM/tstarEM+" m/s")
% % disp("Maneuver 1 state: ["+q1(1)+", "+q1(2)+", "+q1(3)+", "+q1(4)+", "+q1(5)+", "+q1(6)+"]")
% % tauAssist1 = 10*pi;
% % solAssist1 = ode89(odeCR3BPEM, [0 tauAssist1], q1, odeOpts);
% 
% % Deltav2
% q2 = qMan;
% disp("Delta-v 2: "+Deltav2*1000*lstarEM/tstarEM+" m/s")
% disp("Maneuver 2 state: ["+q2(1)+", "+q2(2)+", "+q2(3)+", "+q2(4)+", "+q2(5)+", "+q2(6)+"]")
% tauAssist2 = 2*pi;
% solAssist2 = ode89(odeCR3BPEM, [0 tauAssist2], q2, odeOpts);

% disp("New JC: "+newJC)
% disp("New departure energy:"+newEscE)
% disp("New departure velocity:"+newEscv)

% fig9 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% scatter3(solOrig.y(1,1), solOrig.y(2,1), solOrig.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
% % scatter3(qMan(1), qMan(2), qMan(3), 50, 'w', 'filled', '^', 'DisplayName', "Maneuver")
% % p91 = plot3WithArrows(solOrig.y(1,:), solOrig.y(2,:), solOrig.y(3,:), 'Color', colorMap(flags(idx)+1,:));
% % set(p91, 'DisplayName', "Original Traj.")
% plot3(solOrig.y(1,:), solOrig.y(2,:), solOrig.y(3,:), 'Color', colorMap(flags(idx)+1,:), 'DisplayName', "Original Traj.")
% % % p92 = plot3WithArrows(solAssist1.y(1,:), solAssist1.y(2,:), solAssist1.y(3,:), 'b');
% % % set(p92, 'DisplayName', "Assisted Traj.")
% % plot3(solAssist1.y(1,:), solAssist1.y(2,:), solAssist1.y(3,:), 'b', 'DisplayName', "Assisted Traj.")
% % p93 = plot3WithArrows(solAssist2.y(1,:), solAssist2.y(2,:), solAssist2.y(3,:), 'b');
% % set(p93, 'DisplayName', "Assisted Traj.")
% plot3(solAssist2.y(1,:), solAssist2.y(2,:), solAssist2.y(3,:), 'b', 'DisplayName', "Assisted Traj.")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% % axis([-3 3 -3 3])
% grid on
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg9 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg9.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax9 = gca;
% ax9.SortMethod = 'childorder';
% % exportgraphics(fig9, 'EscapeCR3BP_9.png','BackgroundColor', 'k')

%% Import Cluster Data
% clusterData = load('../PhDScripts/Output/ClusterTrajectoriesCR3BP.mat');
% clusters = clusterData.clusters;
% k = max(clusters);

%% Clustering
% clusterColor = viridis(k);
% cluster = 1;
% 
% fig11 = figure("Position", [200 100 1200 750]);
% hold on
% % for j = 1:nFilt0
% %     q = q0s_filt0(:,j);
% %     tau = tfs_filt0(j);
% %     [~, rout] = mexCR3BP(q, [0 tau], muEM, 1E-12, 1E-12, 1E-10);
% %     if clusters(j) == cluster
% %         plot3(rout(:,1), rout(:,2), rout(:,3), 'Color', [clusterColor(cluster,:), 0.1], 'LineWidth', 1, 'HandleVisibility', 'off')
% %     % else
% %         % plot3(rout(:,1), rout(:,2), rout(:,3), 'Color', [0.5, 0.5, 0.5, 0.05], 'LineWidth', 1, 'HandleVisibility', 'off')
% %     end
% % end
% scatter(q0s_filt0(1,:), q0s_filt0(2,:), 1.75, clusterColor(clusters,:), 'filled', 'HandleVisibility', 'off')
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
% colormap(clusterColor)
% cb11 = colorbar;
% clim([0.5 double(k)+0.5])
% cb11.Ticks = 1:k;
% ylabel(cb11, "Clusters", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb11.Label.Position = cb11.Label.Position+[-1.5 double(k)/2 0];
% leg11 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% leg11.Position = leg11.Position+[0.08 0 0 0];
% drawnow;
% set(leg11.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax11 = gca;
% ax11.SortMethod = 'childorder';
% % exportgraphics(fig11, 'EscapeCR3BP_11.png','BackgroundColor', 'k')

%% Escape Analysis Figure
% nFilt = nFilt1;
% qs = q0s_filt1;
% Es = Es_filt1;
% tfs = tfs_filt1;

% colors = nebula(1000);
% % pointColors = zeros(nFilt0, 3);
% % parfor j = 1:nFilt0
% %     pointColors(j,:) = getColor(colors, ds_filt0(j), [0, ceil(d_max/10000)*10000]);
% % end
% pointColors = zeros(nFilt, 3);
% parfor j = 1:nFilt
%     pointColors(j,:) = getColor(colors, Es(j), [E_min, E_max]);
% end
% % pointColors2 = zeros(nFilt2, 3);
% % parfor j = 1:nFilt2
% %     pointColors2(j,:) = getColor(colors, Es_filt2(j), [E_min, E_max]);
% % end

% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% % for j = 1:nFilt
% %     % if clusters(j) == cluster
% %         q = qs(:,j);
% %         tau = tfs(j);
% %         [~, rout] = mexCR3BP(q, [0 tau], muEM, 1E-12, 1E-12, 1E-10);
% %         % plot3(rout(:,1), rout(:,2), rout(:,3), 'Color', [0.5, 0.5, 0.5, 0.1], 'LineWidth', 1, 'HandleVisibility', 'off')
% %         plot3(rout(:,1), rout(:,2), rout(:,3), 'Color', [pointColors(j,:), 0.1], 'LineWidth', 1, 'HandleVisibility', 'off')
% %         % scatter3(rout(1,1), rout(1,2), rout(1,3), 1.75, clusterColor(cluster,:), 'filled', 'HandleVisibility', 'off')
% %     % end
% % end
% % scatter3(q0s_filt0(1,:), q0s_filt0(2,:), q0s_filt0(3,:), 1.75, color, 'filled', 'HandleVisibility', 'off')
% scatter3(qs(1,:), qs(2,:), qs(3,:), 1.75, pointColors, 'filled', 'HandleVisibility', 'off')
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
% ylabel(cb5, "$\mathcal{E}_{esc}$ [km$^{2}$/s$^{2}$]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb5.Label.Position = cb5.Label.Position+[-3.5 0.135 0];
% % clim([0 ceil(d_max/10000)])
% % cb5.Ticks = 0:1:ceil(d_max/10000);
% % ylabel(cb5, "Perilune [km x$10^{4}$]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% % cb5.Label.Position = cb5.Label.Position+[-2 4.1 0];
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
% volumeFile = 'E:/ApseMapData/CR3BPJCVolume_1_peri_pro_500_2.9_3.17.mat';
% volumeDataFile = 'CR3BPJCVolume_1_peri_pro_500_2.9_3.17.mat';

% volumeFields = who('-file', volumeFile);
% nVolume = length(volumeFields);
% sampleData = load(volumeFile, volumeFields{1});
% sampleFields = fieldnames(sampleData);
% sampleMap = sampleData.(sampleFields{1});
% rSample = size(sampleMap.flags, 1);
% cSample = size(sampleMap.flags, 2);
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
% countsVolume = zeros(nVolume, rSample, cSample);
% flagsVolume = zeros(nVolume, rSample, cSample);
% qVolume = zeros(6, rSample*cSample, nVolume);
% parfor j = 1:nVolume
%     volumeData = load(volumeFile, sortedFields{j});
%     volumeDataFields = fieldnames(volumeData);
%     volumeMap = volumeData.(volumeDataFields{1});
%     countsVolume(j,:) = reshape(volumeMap.counts, 1, rSample, cSample);
%     flagsVolume(j,:) = reshape(volumeMap.flags, 1, rSample, cSample);
%     qVolume(:,:,j) = volumeMap.q;
% end

%% Save JC Volume Data
% save(volumeDataFile, "JCVolume", "countsVolume", "flagsVolume", "qVolume", "-v7.3");

%% Load JC Volume Data
% load(volumeDataFile, "JCVolume", "flagsVolume", "qVolume");
% nVolume = length(JCVolume);
% nSample = size(qVolume, 2);
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
% hist3 = scatter(NaN(1, nSample), NaN(1, nSample), scatSize, repmat(colorMap(7,:), nSample, 1), 'filled', 'HandleVisibility', 'off');
% 
% % xlim(ax3, [-1.25 1.25])
% % ylim(ax3, [-1.25 1.25])
% xlim(ax3, [1-muEM-0.2 1-muEM+0.2])
% ylim(ax3, [-0.2 0.2])
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
% % set(leg3.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% 
% tit3 = annotation(fig3, 'textbox', [0.4 0.93 0.2 0.06], 'String', sprintf('JC = %.4f', JCVolume(1)), 'FontName', 'Times New Roman', 'FontSize', 18, 'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% 
% % v = VideoWriter('EscapeCR3BP_JCVolume_Moon_peri_pro.mp4', 'MPEG-4');
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
% assistDeltav1s = assistedData.Deltav1s;
% assistDeltav2s = assistedData.Deltav2s;
% assistJCs = assistedData.newJCs;
% totalDeltavs = assistDeltav1s+assistDeltav2s;

%% Assisted Escape Figure
% pointColors = zeros(length(totalDeltavs), 3);
% DeltavColors = rdbuInvert(1000);
% % DeltavColors = rdbu(1000);
% maxDeltav = 30;
% % maxDeltav = 50;
% parfor j = 1:n
%     if isnan(totalDeltavs(j))
%         pointColors(j,:) = [1, 1, 1];
%         % pointColors(j,:) = [0, 0, 0];
%     else
%         pointColors(j,:) = getColor(DeltavColors, totalDeltavs(j)*1000*lstarEM/tstarEM, [-maxDeltav, maxDeltav]);        
%     end
% end
% % DeltavColors = copper(1000);
% % maxDeltav = 250;
% % parfor j = 1:n^2
% %     if isnan(maneuverMetrics(j))
% %         pointColors(j,:) = [1, 1, 1];
% %     else
% %         pointColors(j,:) = getColor(DeltavColors, maneuverMetrics(j)*1000, [0, maxDeltav]);        
% %     end
% % end
% 
% fig12 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(xGrid, yGrid, 1.75, pointColors, 'filled', 'HandleVisibility', 'off')
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter(nan, nan, 1.75, 'w', 'filled', 'DisplayName', "Infeasible")
% % scatter(xPeri, yPeri, 30, 'm', 'filled', 'DisplayName', "Manifold Peris.")
% % scatter(sol1.y(1,1), sol1.y(2,1), 50, 'c', 'filled', 'HandleVisibility', 'off')
% scatter(sol2.y(1,1), sol2.y(2,1), 50, 'c', 'filled', 'HandleVisibility', 'off')
% % scatter(trajPeriStates(1,:), trajPeriStates(2,:), 50, 'c', 'o', 'filled', 'HandleVisibility', 'off')
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
% xlabel("$x$ [E-M ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [E-M ndim]", 'Interpreter', 'latex')
% title("Earth-Moon Rot.: JC = "+JC, 'Interpreter', 'latex')
% colormap(rdbuInvert)
% % colormap(rdbu)
% cb12 = colorbar;
% clim([-maxDeltav maxDeltav])
% cb12.Ticks = -maxDeltav:10:maxDeltav;
% cb12.TickLabels = [string(-maxDeltav)+"+", string(-(maxDeltav-10):10:maxDeltav-10), string(maxDeltav)+"+"];
% ylabel(cb12, "$\Delta v$ [m/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb12.Label.Position = cb12.Label.Position+[-3 31 0];
% % cb12.Label.Position = cb12.Label.Position+[-3 51 0];
% % colormap(copper)
% % cb12 = colorbar;
% % clim([0 maxDeltav])
% % cb12.Ticks = 0:50:maxDeltav;
% % cb12.TickLabels = [0:50:maxDeltav-50, string(maxDeltav)+"+"];
% % ylabel(cb12, "$\Delta v$ Saved [m/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% % cb12.Label.Position = cb12.Label.Position+[-3 130 0];
% leg12 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% leg12.Position = leg12.Position+[0.1 0 0 0];
% drawnow;
% set(leg12.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% ax12 = gca;
% ax12.SortMethod = 'childorder';
% % exportgraphics(fig12, 'EscapeCR3BP_12.png', 'BackgroundColor', 'k')
% % exportgraphics(fig12, 'EscapeCR3BP_12.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Test Trajectory
% xSample = 0.0125251;
% ySample = -0.473447;
% idx = find((abs(xGrid-xSample) < 1E-5) & (abs(yGrid-ySample) < 1E-5));
% qOrig = map.q(:,idx);
% disp("Original IC: ["+qOrig(1)+", "+qOrig(2)+", "+qOrig(3)+", "+qOrig(4)+", "+qOrig(5)+", "+qOrig(6)+"]")
% tauOrig = 4*pi;
% solOrig = ode89(odeCR3BPEM, [0 tauOrig], qOrig, odeOpts);
% vOrig = norm(qOrig(4:5));
% vhat = qOrig(4:5)./vOrig;
% qAssist = qOrig;
% Deltav = Deltav1s(idx);
% qAssist(4:5) = (vOrig+Deltav).*vhat;
% disp("Delta-v: "+Deltav*1000*lstarEM/tstarEM+" m/s")
% disp("Assisted IC: ["+qAssist(1)+", "+qAssist(2)+", "+qAssist(3)+", "+qAssist(4)+", "+qAssist(5)+", "+qAssist(6)+"]")
% tauAssist = 4*pi;
% solAssist = ode89(odeCR3BPEM, [0 tauAssist], qAssist, odeOpts);

% fig7 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% scatter3(solOrig.y(1,1), solOrig.y(2,1), solOrig.y(3,1), 50, 'g', 'filled', 'DisplayName', "Start")
% % p71 = plot3WithArrows(solOrig.y(1,:), solOrig.y(2,:), solOrig.y(3,:), 'Color', colorMap(flags(idx)+1,:));
% % set(p71, 'DisplayName', "Original Traj.")
% plot3(solOrig.y(1,:), solOrig.y(2,:), solOrig.y(3,:), 'Color', colorMap(flags(idx)+1,:), 'DisplayName', "Original Traj.")
% % p72 = plot3WithArrows(solAssist.y(1,:), solAssist.y(2,:), solAssist.y(3,:), 'g');
% % set(p72, 'DisplayName', "Assisted Traj.")
% plot3(solAssist.y(1,:), solAssist.y(2,:), solAssist.y(3,:), 'g', 'DisplayName', "Assisted Traj.")
% axis equal
% axis([-1.25 1.25 -1.25 1.25])
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