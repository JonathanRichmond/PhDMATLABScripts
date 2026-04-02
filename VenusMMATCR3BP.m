%%% VenusMMATCR3BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 31 February 2026
%%% U: 1 April 2026

clear

%NEED TO FIX TRANSFER ANGLES
%% Import Transfer Data
transferData = load('E:/MMATData/VenusMMATCR3BP_2_97_L2Butterfly_fliptrue_3_0003_L2Halo_flipfalse.mat');
transfers = cellfun(@(n) transferData.(n), fieldnames(transferData));
n_transfers = length(transfers);
disp("Total transfers: "+n_transfers)
TOFs = zeros(n_transfers,1);
Deltav_1s = zeros(n_transfers,1);
Deltav_2s = zeros(n_transfers,1);
theta_E_deps = zeros(n_transfers,1);
for j = 1:length(transfers)
    currentTransfer = transfers(j);
    TOFs(j) = currentTransfer.TOF;
    Deltav_1s(j) = currentTransfer.Deltav_1;
    Deltav_2s(j) = currentTransfer.Deltav_2;
    theta_E_deps(j) = wrapTo360(atan2d(currentTransfer.departureConic.state(2), currentTransfer.departureConic.state(1)));
end
Deltavs = Deltav_1s+Deltav_2s;
[~, sortDeltavs]= sort(Deltavs);
[~, sortTOFs] = sort(TOFs);
disp("Min. Deltav: "+Deltavs(sortDeltavs(1))+" km/s")
disp("Min. TOF: "+TOFs(sortTOFs(1))/24/3600/365.25+" yrs")
transferIdx = sortDeltavs(30);
transfer = transfers(transferIdx);

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

%% Sun-Venus Data
gmV = 3.2485859200000000E5; % Venus gravitational parameter [km^3/s^2]
mV = gmV/6.67384E-20; % Venus mass [kg]
RV = 6.0518000000000002E3; % Venus radius [km]

muSV = gmV/(gmS+gmV); % Mass ratio
mstarSV = (gmS+gmV)/6.67384E-20; % Characteristic mass [kg]
lstarSV = 1.0820891506628585E8; % Characteristic length [km]
tstarSV = sqrt(lstarSV^3/(gmS+gmV)); % Characteristic time [s]

g1SV = muSV; % Initial guess
delg1 = 1;
while abs(delg1) > eps
    f = ((1-muSV)/((1-g1SV)^2))-(muSV/(g1SV^2))-1+muSV+g1SV;
    fprime = ((2*(1-muSV))/((1-g1SV)^3))+((2*muSV)/(g1SV^3))+1;
    g1SVnew = g1SV-(f/fprime);
    delg1 = g1SVnew-g1SV;
    g1SV = g1SVnew;
    i = i+1;
end
a1SV = 1-muSV-g1SV;

g2SV = muSV; % Initial guess
delg2 = 1;
while abs(delg2) > eps
    f = ((1-muSV)/((1+g2SV)^2))+(muSV/(g2SV^2))-1+muSV-g2SV;
    fprime = ((-2*(1-muSV))/((1+g2SV)^3))-((2*muSV)/(g2SV^3))-1;
    g2SVnew = g2SV-(f/fprime);
    delg2 = g2SVnew-g2SV;
    g2SV = g2SVnew;
    i = i+1;
end
a2SV = 1-muSV+g2SV;

g3SV = muSV; %Initial guess
delg3 = 1;
while abs(delg3) > eps
    f = ((1-muSV)/(g3SV^2))+(muSV/((1+g3SV)^2))-muSV-g3SV;
    fprime = ((-2*(1-muSV))/(g3SV^3))-((2*muSV)/((1+g3SV)^3))-1;
    g3SVnew = g3SV-(f/fprime);
    delg3 = g3SVnew-g3SV;
    g3SV = g3SVnew;
    i = i+1;
end
a3SV = -1*muSV-g3SV;

a45SV = 0.5-muSV;
b4SV = sqrt(3)/2;
b5SV = -b4SV;

%% Design Variables
RSoIm = lstarSE*(mm/mS)^(2/5); % Moon sphere of influence radius [km]
RSoIE = 0.09877*lstarSE; % Earth sphere of influence radius [km]
RSoIV = 0.1107*lstarSV; % Venus sphere of influence radius [km]

%% Frame Transformations
e_dep_trans = transfer.t_0+transfer.departureManifoldArc1.TOF*tstarEM;
e_dep_SoI = e_dep_trans+transfer.departureManifoldArc2.TOF*tstarSE;
e_bridge = e_dep_SoI+transfer.departureConic.TOF;
e_arr_int = e_bridge+transfer.bridgeConic.TOF;

theta_rot = deg2rad(theta_E_deps(transferIdx));
R = [cos(theta_rot), sin(theta_rot) 0; -sin(theta_rot), cos(theta_rot) 0; 0 0 1];

t_EM = transfer.t_0/tstarEM+transfer.departureManifoldArc1.t;
q_EM = [transfer.departureManifoldArc1.x, transfer.departureManifoldArc1.y, transfer.departureManifoldArc1.z, transfer.departureManifoldArc1.xdot, transfer.departureManifoldArc1.ydot, transfer.departureManifoldArc1.zdot];
q_EI_EM = rotToP1EclipJ2000(muEM, transfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, t_EM, q_EM);
q_EI_SE = [q_EI_EM(:,1:3).*lstarEM./lstarSE, q_EI_EM(:,4:6).*lstarEM.*tstarSE./lstarSE./tstarEM];
t_EM_SE = t_EM.*tstarEM./tstarSE;
q_EM_SE = P2EclipJ2000ToRot(muSE, transfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, t_EM_SE, q_EI_SE);
q_EM_SI = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_EM_SE, q_EM_SE);
Q_EM_SI = [q_EM_SI(:,1:3).*lstarSE, q_EM_SI(:,4:6).*lstarSE./tstarSE];
r_EM_SI_R = (R*Q_EM_SI(:,1:3)')';

t_SE = e_dep_trans/tstarSE+transfer.departureManifoldArc2.t;
q_SE = [transfer.departureManifoldArc2.x, transfer.departureManifoldArc2.y, transfer.departureManifoldArc2.z, transfer.departureManifoldArc2.xdot, transfer.departureManifoldArc2.ydot, transfer.departureManifoldArc2.zdot];
q_SE_SI = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_SE, q_SE);
Q_SE_SI = [q_SE_SI(:,1:3).*lstarSE, q_SE_SI(:,4:6).*lstarSE./tstarSE];
r_SE_SI_R = (R*Q_SE_SI(:,1:3)')';

ode = @(t,r) ODE_C(t, r, gmS);
odeOpts = odeset('RelTol', 1E-12, 'AbsTol', 1E-12);
tspan_dep = e_dep_SoI+[0:(3600*24):transfer.departureConic.TOF, transfer.departureConic.TOF];
IC_dep = transfer.departureConic.state;
dep_sol = ode89(ode, tspan_dep, IC_dep, odeOpts);
r_dep_R = R*dep_sol.y(1:3,:);

tspan_bridge = e_bridge+[0:(3600*24):transfer.bridgeConic.TOF, transfer.bridgeConic.TOF];
IC_bridge = transfer.bridgeConic.state;
bridge_sol = ode89(ode, tspan_bridge, IC_bridge, odeOpts);
r_bridge_R = R*bridge_sol.y(1:3,:);

tspan_arr = e_arr_int+[0:(3600*24):transfer.arrivalConic.TOF, transfer.arrivalConic.TOF];
IC_arr = transfer.arrivalConic.state;
arr_sol = ode89(ode, tspan_arr, IC_arr, odeOpts);
r_arr_R = R*arr_sol.y(1:3,:);

t_SV = transfer.theta_arr_f+transfer.arrivalManifoldArc.t; % Because Venus is not yet realistically phased
q_SV = [transfer.arrivalManifoldArc.x, transfer.arrivalManifoldArc.y, transfer.arrivalManifoldArc.z, transfer.arrivalManifoldArc.xdot, transfer.arrivalManifoldArc.ydot, transfer.arrivalManifoldArc.zdot];
q_SV_SI = planetRotToSunEclipJ2000(muSV, transfer.initialEpoch, 'Venus', lstarSV, tstarSV, t_SV, q_SV);
Q_SV_SI = [q_SV_SI(:,1:3).*lstarSV, q_SV_SI(:,4:6).*lstarSV./tstarSV];
r_SV_SI_R = (R*Q_SV_SI(:,1:3)')';

t_m_EM = linspace(transfer.t_0/tstarEM, e_dep_SoI/tstarEM, 10001);
q_m_EI_EM = rotToP1EclipJ2000(muEM, transfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, t_m_EM, ones(10001, 1)*[1-muEM, 0, 0, 0, 0, 0]);
q_m_EI_SE = [q_m_EI_EM(:,1:3).*lstarEM./lstarSE, q_m_EI_EM(:,4:6).*lstarEM.*tstarSE./lstarSE./tstarEM];
t_m_SE = t_m_EM.*tstarEM./tstarSE;
q_m_SE = P2EclipJ2000ToRot(muSE, transfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, t_m_SE, q_m_EI_SE);

t_E_SE = linspace(transfer.t_0/tstarSE, transfer.t_0/tstarSE+2*pi, 10001);
q_E_SI = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_E_SE, ones(10001, 1)*[1-muSE, 0, 0, 0, 0, 0]);
Q_E_SI = [q_E_SI(:,1:3).*lstarSE, q_E_SI(:,4:6).*lstarSE./tstarSE];

t_V_SV = linspace(transfer.theta_arr_f-2*pi, transfer.theta_arr_f, 10001);
q_V_SI = planetRotToSunEclipJ2000(muSV, transfer.initialEpoch, 'Venus', lstarSV, tstarSV, t_V_SV, ones(10001, 1)*[1-muSV, 0, 0, 0, 0, 0]);
Q_V_SI = [q_V_SI(:,1:3).*lstarSV, q_V_SI(:,4:6).*lstarSV./tstarSV];

%% Earth-Moon Trajectory Plot
% fig1 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM, [1-muEM, 0, 0]);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% p11 = plot3WithArrows(transfer.departureOrbit.x, transfer.departureOrbit.y, transfer.departureOrbit.z, 'g--', 'NumArrows', 2, 'ArrowScale', 2);
% set(p11, 'DisplayName', "Dep. Orbit")
% p12 = plot3WithArrows(transfer.departureManifoldArc1.x, transfer.departureManifoldArc1.y, transfer.departureManifoldArc1.z, 'r');
% set(p12, 'DisplayName', "Dep. CR3BP Arc")
% plot3(1-muEM+RSoIm/lstarEM*cos(linspace(0, 2*pi, 101)), RSoIm/lstarEM*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Moon SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
% zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
% % title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg1.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig1, 'VenusMMATCR3BP_1.png', 'BackgroundColor', 'k')

%% Sun-Earth Trajectory Plot
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", 10*RE/lstarSE, [1-muSE, 0, 0]);
% set(Earth, 'DisplayName', "Earth (x10)")
% p21 = plot3WithArrows(q_m_SE(:,1), q_m_SE(:,2), q_m_SE(:,3), 'w--');
% set(p21, 'DisplayName', "Moon Traj")
% scatter3(a1SE, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "SE $L_{1}$")
% % scatter3(a2SE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "SE $L_{2}$")
% p22 = plot3WithArrows(q_EM_SE(:,1), q_EM_SE(:,2), q_EM_SE(:,3), 'r', 'NumArrows', 2, 'ArrowScale', 2);
% set(p22, 'DisplayName', "Dep. CR3BP Arc")
% p23 = plot3WithArrows(transfer.departureManifoldArc2.x(1:end-7), transfer.departureManifoldArc2.y(1:end-7), transfer.departureManifoldArc2.z(1:end-7), 'r', 'ArrowScale', 0.5);
% set(p23, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% % plot3(1-muSE+RSoIE/lstarSE*cos(linspace(0, 2*pi, 101)), RSoIE/lstarSE*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Earth SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [SE ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [SE ndim]", 'Interpreter', 'latex')
% hz2 = zlabel("$z$ [SE ndim]", 'Interpreter', 'latex');
% % shiftZLabel(gca, hz2)
% % title("Sun-Earth Rot.", 'Interpreter', 'latex')
% leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg2.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig2, 'VenusMMATCR3BP_2.png', 'BackgroundColor', 'k')

%% Heliocentric MMAT Plot
% disp("TOF: "+transfer.TOF/24/3600/365.25+" yrs")

% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% Sun = plot3DBody("Sun", 10*RS, [0, 0, 0]);
% set(Sun, 'DisplayName', "Sun (x10)")
% p31 = plot3WithArrows(Q_E_SI(:,1), Q_E_SI(:,2), Q_E_SI(:,3), 'g--', 'NumArrows', 2, 'ArrowScale', 0.75);
% set(p31, 'DisplayName', "Earth Orbit")
% p32 = plot3WithArrows(Q_V_SI(:,1), Q_V_SI(:,2), Q_V_SI(:,3), '--', 'Color', [1 0.5 0], 'NumArrows', 2);
% set(p32, 'DisplayName', "Venus Orbit")
% Earth = plot3DBody("Earth", 1000*RE, Q_E_SI(1,1:3));
% set(Earth, 'DisplayName', "Earth (x1000) at Dep.")
% Venus = plot3DBody("Venus", 1000*RV, Q_V_SI(end,1:3));
% set(Venus, 'DisplayName', "Venus (x1000) at Arr.")
% p33 = plot3WithArrows(Q_EM_SI(:,1), Q_EM_SI(:,2), Q_EM_SI(:,3), 'r', 'NumArrows', 2);
% set(p33, 'DisplayName', "Dep. CR3BP Arc")
% p34 = plot3WithArrows(Q_SE_SI(:,1), Q_SE_SI(:,2), Q_SE_SI(:,3), 'r');
% set(p34, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% p35 = plot3WithArrows(dep_sol.y(1,:), dep_sol.y(2,:), dep_sol.y(3,:), 'm', 'ArrowScale', 0.75);
% set(p35, 'DisplayName', "Dep. Conic")
% p36 = plot3WithArrows(bridge_sol.y(1,:), bridge_sol.y(2,:), bridge_sol.y(3,:), 'Color', [0.5 0 0.5]);
% set(p36, 'DisplayName', "Bridge Conic")
% p37 = plot3WithArrows(arr_sol.y(1,:), arr_sol.y(2,:), arr_sol.y(3,:), 'c');
% set(p37, 'DisplayName', "Arr. Conic")
% p38 = plot3WithArrows(Q_SV_SI(:,1), Q_SV_SI(:,2), Q_SV_SI(:,3), 'b', 'FlipDir', 'On');
% set(p38, 'DisplayName', "Arr. CR3BP Arc")
% scatter3(bridge_sol.y(1,1), bridge_sol.y(2,1), bridge_sol.y(3,1), 75, 'w', 'filled', 's', 'DisplayName', "$\Delta v_{1}="+num2str(transfer.Deltav_1)+"$ km/s")
% scatter3(bridge_sol.y(1,end), bridge_sol.y(2,end), bridge_sol.y(3,end), 75, 'w', 'filled', '^', 'DisplayName', "$\Delta v_{2}="+num2str(transfer.Deltav_2)+"$ km/s")
% axis equal
% grid on
% xlabel("$X$ [km]", 'Interpreter', 'latex')
% ylabel("$Y$ [km]", 'Interpreter', 'latex')
% hz3 = zlabel("$Z$ [km]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz3)
% % title("Sun-Centered Ecliptic J2000 $|$ TOF = "+num2str(transfer.TOF/24/3600/365.25, 3)+" yrs", 'Interpreter', 'latex')
% leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg3.EntryContainer.NodeChildren(end-3).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig3, 'VenusMMATCR3BP_3.png', 'BackgroundColor', 'k')

% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% Sun = plot3DBody("Sun", 10*RS, [0, 0, 0]);
% set(Sun, 'DisplayName', "Sun (x10)")
% p31 = plot3WithArrows(Q_E_SI(:,1), Q_E_SI(:,2), Q_E_SI(:,3), 'g--', 'NumArrows', 2, 'ArrowScale', 0.75);
% set(p31, 'DisplayName', "Earth Orbit")
% p32 = plot3WithArrows(Q_V_SI(:,1), Q_V_SI(:,2), Q_V_SI(:,3), '--', 'Color', [1 0.5 0], 'NumArrows', 2);
% set(p32, 'DisplayName', "Venus Orbit")
% Earth = plot3DBody("Earth", 1000*RE, [lstarSE, 0, 0]);
% set(Earth, 'DisplayName', "Earth (x1000) at Dep. Conic")
% Venus = plot3DBody("Venus", 1000*RV, (R*Q_V_SI(end,1:3)')');
% set(Venus, 'DisplayName', "Venus (x1000) at Arr.")
% p33 = plot3WithArrows(r_EM_SI_R(:,1), r_EM_SI_R(:,2), r_EM_SI_R(:,3), 'r', 'NumArrows', 2);
% set(p33, 'DisplayName', "Dep. CR3BP Arc")
% p34 = plot3WithArrows(r_SE_SI_R(:,1), r_SE_SI_R(:,2), r_SE_SI_R(:,3), 'r');
% set(p34, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% p35 = plot3WithArrows(r_dep_R(1,:), r_dep_R(2,:), r_dep_R(3,:), 'm', 'ArrowScale', 0.75);
% set(p35, 'DisplayName', "Dep. Conic")
% p36 = plot3WithArrows(r_bridge_R(1,:), r_bridge_R(2,:), r_bridge_R, 'Color', [0.5 0 0.5]);
% set(p36, 'DisplayName', "Bridge Conic")
% p37 = plot3WithArrows(r_arr_R(1,:), r_arr_R(2,:), r_arr_R(3,:), 'c');
% set(p37, 'DisplayName', "Arr. Conic")
% p38 = plot3WithArrows(r_SV_SI_R(:,1), r_SV_SI_R(:,2), r_SV_SI_R(:,3), 'b', 'FlipDir', 'On');
% set(p38, 'DisplayName', "Arr. CR3BP Arc")
% scatter3(r_bridge_R(1,1), r_bridge_R(2,1), r_bridge_R(3,1), 75, 'w', 'filled', 's', 'DisplayName', "$\Delta v_{1}="+num2str(transfer.Deltav_1)+"$ km/s")
% scatter3(r_bridge_R(1,end), r_bridge_R(2,end), r_bridge_R(3,end), 75, 'w', 'filled', '^', 'DisplayName', "$\Delta v_{2}="+num2str(transfer.Deltav_2)+"$ km/s")
% axis equal
% grid on
% xlabel("$X$ [km]", 'Interpreter', 'latex')
% ylabel("$Y$ [km]", 'Interpreter', 'latex')
% % title("Sun-Centered Ecliptic J2000 $|$ TOF = "+num2str(transfer.TOF/24/3600/365.25, 3)+" yrs", 'Interpreter', 'latex')
% leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg3.EntryContainer.NodeChildren(end-3).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig3, 'VenusMMATCR3BP_3.png', 'BackgroundColor', 'k')

%% Sun-Venus Trajectory Plot
% fig4 = figure("Position", [200 100 1200 750]);
% hold on
% Venus = plot3DBody("Venus", 10*RV/lstarSV, [1-muSV, 0, 0]);
% set(Venus, 'DisplayName', "Venus (x10)")
% % scatter3(a1SV, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "SV $L_{1}$")
% scatter3(a2SV, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "SV $L_{2}$")
% p41 = plot3WithArrows(transfer.arrivalOrbit.x, transfer.arrivalOrbit.y, transfer.arrivalOrbit.z, 'g--', 'NumArrows', 2, 'ArrowScale', 1.5);
% set(p41, 'DisplayName', "Arr. Orbit")
% p42 = plot3WithArrows(transfer.arrivalManifoldArc.x(1:end-6), transfer.arrivalManifoldArc.y(1:end-6), transfer.arrivalManifoldArc.z(1:end-6), 'b', 'FlipDir', 'On', 'ArrowScale', 0.5);
% set(p42, 'DisplayName', "Arr. CR3BP Arc")
% % plot3(1-muSV+RSoIV/lstarSV*cos(linspace(0, 2*pi, 101)), RSoIV/lstarSV*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Venus SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [SV ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [SV ndim]", 'Interpreter', 'latex')
% hz4 = zlabel("$z$ [SV ndim]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz4)
% % title("Sun-Venus Rot.", 'Interpreter', 'latex')
% leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig4, 'VenusMMATCR3BP_4.png', 'BackgroundColor', 'k')

%% Family Plots
years = 0; % 2030
DeltavHohmann = 5.954; % km/s (Venus)

% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(wrapTo360(theta_E_0s), TOFs/24/3600/365.25, 20, Deltavs, 'filled', 'HandleVisibility', 'off')
% xlim([0 360])
% grid on
% xlabel("$\theta_{E,0}$ [deg]", 'Interpreter', 'latex')
% ylabel("TOF [yrs]", 'Interpreter', 'latex')
% % title("MMAT Family Tradespace", 'Interpreter', 'latex')
% colormap(viridis)
% cb5 = colorbar;
% caxis([4 10])
% ylabel(cb5, "$\Delta v$ [km/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb5.Label.Position = cb5.Label.Position+[-2 3.1 0];
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig5, 'VenusMMATCR3BP_5.png', 'BackgroundColor', 'k')

% fig6 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(wrapTo360(((6.253075709008801+years*2*pi)*180/pi):0.5:((6.253075709008801+(years+1)*2*pi)*180/pi)), wrapTo360(((5.679677667107661+years*2*pi*tstarSE/tstarSV)*180/pi):(0.5*tstarSE/tstarSV):((5.679677667107661+(years+1)*2*pi*tstarSE/tstarSV)*180/pi)), 5, 'w', 'filled', 'HandleVisibility', 'off')
% scatter(wrapTo360(theta_E_0s), wrapTo360(theta_V_0s), 20, Deltavs, 'filled', 'HandleVisibility', 'off')
% axis equal
% axis([0 360 0 360])
% grid on
% xlabel("$\theta_{E,0}$ [deg]", 'Interpreter', 'latex')
% ylabel("$\theta_{V,0}$ [deg]", 'Interpreter', 'latex')
% % title("MMAT Family Tradespace", 'Interpreter', 'latex')
% colormap(viridis)
% cb6 = colorbar;
% caxis([4 10])
% ylabel(cb6, "$\Delta v$ [km/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb6.Label.Position = cb6.Label.Position+[-2 3.1 0];
% set(gca, 'Color', 'k');
% view(2)
% hold off
% % exportgraphics(fig6, 'VenusMMATCR3BP_6.png', 'BackgroundColor', 'k')

%% Comparison Metric
mask = Deltavs < DeltavHohmann;
DeltavsComp = Deltavs(mask);
TOFsComp = TOFs(mask);
[sortTOFsComp, sortTOFsCompIdx] = sort(TOFsComp);
sortDeltavsComp = DeltavsComp(sortTOFsCompIdx);

sortParetoMask = false(size(sortTOFsComp));
bestDeltav = inf;
for j = 1:length(sortTOFsComp)
    if sortDeltavsComp(j) < bestDeltav
        sortParetoMask(j) = true;
        bestDeltav = sortDeltavsComp(j);
    end
end
paretoMask = false(size(TOFsComp));
paretoMask(sortTOFsCompIdx) = sortParetoMask;
TOFsPareto = TOFsComp(paretoMask);
DeltavsPareto = DeltavsComp(paretoMask);
[sortTOFsPareto, sortTOFsParetoIdx] = sort(TOFsPareto);
sortDeltavsPareto = DeltavsPareto(sortTOFsParetoIdx);

pp = pchip(sortTOFsPareto/24/3600/365.25, sortDeltavsPareto);
TOFFine = linspace(min(sortTOFsPareto)/24/3600/365.25, max(sortTOFsPareto)/24/3600/365.25, 500);
DeltavSmooth = ppval(pp, TOFFine);

refTOF = 5;
TOFMaxPareto = max(sortTOFsPareto)/24/3600/365.25;
TOFGrid = linspace(sortTOFsPareto(1)/24/3600/365.25, refTOF, 1000);
DeltavGrid = ppval(pp, TOFGrid);
DeltavGrid(TOFGrid > TOFMaxPareto) = sortDeltavsPareto(end);
HV = trapz(TOFGrid, DeltavHohmann*ones(size(DeltavGrid))-DeltavGrid);
disp("Hypervolume: "+HV+" yrs*km/s")

DeltavRef = DeltavHohmann*ones(length(TOFsPareto)+1,1);
xFill = [[sortTOFsPareto/24/3600/365.25; refTOF]; flipud([sortTOFsPareto/24/3600/365.25; refTOF])];
yFill = [DeltavRef; flipud([sortDeltavsPareto; sortDeltavsPareto(end)])];

TOFs_n = (TOFsComp-min(TOFsComp))/(max(TOFsComp)-min(TOFsComp));
Deltavs_n = (DeltavsComp-min(DeltavsComp))/(max(DeltavsComp)-min(DeltavsComp));
TOFsPareto_n = (TOFsPareto-min(TOFsPareto))/(max(TOFsPareto)-min(TOFsPareto));
DeltavsPareto_n = (DeltavsPareto-min(DeltavsPareto))/(max(DeltavsPareto)-min(DeltavsPareto));
pareto = [TOFsPareto_n, DeltavsPareto_n];
points = [TOFs_n, Deltavs_n];
distPareto = zeros(size(points, 1),1);
for j = 1:size(points, 1)
    d = sqrt(sum((pareto-points(j,:)).^2, 2));
    distPareto(j) = min(d);
end
DM = mean(distPareto);
disp("Mean distance: "+DM)

xgrid = linspace(min(TOFsComp)/24/3600/365.25, max(TOFsComp)/24/3600/365.25, 100);
ygrid = linspace(min(DeltavsComp), max(DeltavsComp), 100);
[Xg, Yg] = meshgrid(xgrid, ygrid);
data = [TOFsComp(:)/24/3600/365.25, DeltavsComp(:)];
gridPoints = [Xg(:), Yg(:)];
[f, ~] = ksdensity(data, gridPoints);
F = reshape(f, size(Xg));

fig7 = figure("Position", [200 100 1200 750]);
hold on
scatter(TOFs/24/3600/365.25, Deltavs, 10, 'filled', 'MarkerFaceAlpha', 0.25, 'DisplayName', "Transfers")
% scatter(TOFs/24/3600/365.25, Deltavs, 10, angleColor(deg2rad(transferAnglesA)), 'filled', 'HandleVisibility', 'off')
% scatter(TOFs/24/3600/365.25, Deltavs, 10, angleColor(deg2rad(transferAnglesB)), 'filled', 'HandleVisibility', 'off')
% scatter(TOFs/24/3600/365.25, Deltavs, 10, angleColor(deg2rad(transferAnglesA-transferAnglesB)), 'filled', 'HandleVisibility', 'off')
% scatter(TOFs(transferIdx)/24/3600/365.25, Deltavs(transferIdx), 30, 'w', 'filled', 'HandleVisibility', 'off')
yline(DeltavHohmann, 'k--', 'DisplayName', "Hohmann Baseline")
plot(TOFFine, DeltavSmooth, 'Color', [0.75 0 0], 'DisplayName', "Pareto Front")
fill(xFill, yFill, [1 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', "Hypervolume")
contour(Xg, Yg, F, 5, 'LineColor', 0.25*[1 1 1], 'LineWidth', 1, 'DisplayName', 'Density Contours')
axis([1 5 3 10])
grid on
xlabel("TOF [yrs]", 'Interpreter', 'latex')
ylabel("$\Delta v$ [km/s]", 'Interpreter', 'latex')
% title("$3.07$ Southern $L_{2}$ Halo - $3.0003$ Northern $L_{2}$ Halo", 'Interpreter', 'latex')
% leg7 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% phasemap;
% pb7 = phasebar('deg', 'Location', 'northeast', 'Size', 0.275, 'Title', "Test");
% tpb7 = title(pb7, '$\phi$ [deg]', 'Interpreter', 'latex');
% tpb7 = title(pb7, '$\phi_{H}$ [deg]', 'Interpreter', 'latex');
% tpb7 = title(pb7, '$\phi_{M}$ [deg]', 'Interpreter', 'latex');
% tpb7.Position(2) = tpb7.Position(2)-40;
set(gca, 'Color', 'w');
view(2)
hold off
% exportgraphics(fig7, 'VenusMMATCR3BP_7.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Family Comparison
% [142:167, 169:194, 196:221]
sheetData = readcell('Transfer Metrics.xlsx', 'Sheet', 'Venus Orientation', 'Range', 'A196:L221');
JC_Lyap = cell2mat(sheetData(6:11,3));
JC_halo = cell2mat(sheetData(1:4,3));
JC_vert = cell2mat(sheetData(13:18,3));
JC_axial = cell2mat(sheetData(20:21,3));
JC_butt = cell2mat(sheetData(23:26,3));
HV_Lyap = cell2mat(sheetData(6:11,11));
HV_halo = cell2mat(sheetData(1:4,11));
HV_vert = cell2mat(sheetData(13:18,11));
HV_axial = cell2mat(sheetData(20:21,11));
HV_butt = cell2mat(sheetData(23:26,11));
MD_Lyap = cell2mat(sheetData(6:11,12));
MD_halo = cell2mat(sheetData(1:4,12));
MD_vert = cell2mat(sheetData(13:18,12));
MD_axial = cell2mat(sheetData(20:21,12));
MD_butt = cell2mat(sheetData(23:26,12));

fig9 = figure("Position", [200 100 1200 750]);
tiledlayout(2,1)
ax1 = nexttile;
hold on
plot(JC_Lyap, HV_Lyap, 'r*-', 'DisplayName', "$L_{2}$ Lyapunov")
plot(JC_halo, HV_halo, 'g*-', 'DisplayName', "$L_{2}$ N Halo")
plot(JC_vert, HV_vert, 'b*-', 'DisplayName', "$L_{2}$ Vertical")
plot(JC_axial, HV_axial, 'm*-', 'DisplayName', "$L_{2}$ SW Axial")
plot(JC_butt, HV_butt, 'c*-', 'DisplayName', "$L_{2}$ N Butterfly")
xticks([2.97, 3.0, 3.03, 3.07, 3.1, 3.13])
grid on
ylabel("Hypervolume", 'Interpreter', 'latex')
leg9 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
hold off
ax2 = nexttile;
hold on
plot(JC_Lyap, MD_Lyap, 'r*-', 'DisplayName', "$L_{2}$ Lyapunov")
plot(JC_halo, MD_halo, 'g*-', 'DisplayName', "$L_{2}$ N Halo")
plot(JC_vert, MD_vert, 'b*-', 'DisplayName', "$L_{2}$ Vertical")
plot(JC_axial, MD_axial, 'm*-', 'DisplayName', "$L_{2}$ SW Axial")
plot(JC_butt, MD_butt, 'c*-', 'DisplayName', "$L_{2}$ N Butterfly")
xticks([2.97, 3.0, 3.03, 3.07, 3.1, 3.13])
grid on
ylabel("Mean Distance", 'Interpreter', 'latex')
hold off
xlabel("Earth-Moon JC", 'Interpreter', 'latex')
% title("Earth-Moon Departure Orbit - Sun-Venus $3.000713$ $L_{2}$ Southern Halo", 'Interpreter', 'latex')
set(gca, 'Color', 'w');
% exportgraphics(fig9, 'VenusMMATCR3BP_9.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')