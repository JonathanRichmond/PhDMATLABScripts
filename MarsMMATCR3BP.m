%%% MarsMMATCR3BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 31 January 2026
%%% U: 18 February 2026

clear

transferData = load('../PhDScripts/Output/MMAT/MarsMMATCR3BP_3_03_L1Halo_3_0001857_L1Halo.mat');
transfers = cellfun(@(n) transferData.(n), fieldnames(transferData));
n_transfers = length(transfers);
disp("Total transfers: "+n_transfers)
TOFs = zeros(n_transfers,1);
Deltav_1s = zeros(n_transfers,1);
Deltav_2s = zeros(n_transfers,1);
theta_E_0s = zeros(n_transfers,1);
theta_M_0s = zeros(n_transfers,1);
for j = 1:length(transfers)
    TOFs(j) = transfers(j).TOF;
    Deltav_1s(j) = transfers(j).Deltav_1;
    Deltav_2s(j) = transfers(j).Deltav_2;
    theta_E_0s(j) = transfers(j).theta_dep_0*180/pi;
    theta_M_0s(j) = transfers(j).theta_arr_0*180/pi;
end
Deltavs = Deltav_1s+Deltav_2s;
[~, sortDeltavs]= sort(Deltavs);
[~, sortTOFs] = sort(TOFs);
disp("Min. Deltav: "+Deltavs(sortDeltavs(1))+" km/s")
disp("Min. TOF: "+TOFs(sortTOFs(1))/24/3600/365.25+" yrs")
transfer = transfers(1);

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

%% Design Variables
RSoIm = lstarSE*(mm/mS)^(2/5); % Moon sphere of influence radius [km]
RSoIE = 0.09877*lstarSE; % Earth sphere of influence radius [km]
RSoIM = 0.05375*lstarSM; % Mars sphere of influence radius [km]

%% Frame Transformations
e_dep_trans = transfer.t_0+transfer.departureManifoldArc1.TOF*tstarEM;
e_dep_SoI = e_dep_trans+transfer.departureManifoldArc2.TOF*tstarSE;
e_bridge = e_dep_SoI+transfer.departureConic.TOF;
e_arr_int = e_bridge+transfer.bridgeConic.TOF;

t_EM = transfer.t_0/tstarEM+transfer.departureManifoldArc1.t;
q_EM = [transfer.departureManifoldArc1.x, transfer.departureManifoldArc1.y, transfer.departureManifoldArc1.z, transfer.departureManifoldArc1.xdot, transfer.departureManifoldArc1.ydot, transfer.departureManifoldArc1.zdot];
q_EI_EM = rotToP1EclipJ2000(muEM, transfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, t_EM, q_EM);
q_EI_SE = [q_EI_EM(:,1:3).*lstarEM./lstarSE, q_EI_EM(:,4:6).*lstarEM.*tstarSE./lstarSE./tstarEM];
t_EM_SE = t_EM.*tstarEM./tstarSE;
q_EM_SE = P2EclipJ2000ToRot(muSE, transfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, t_EM_SE, q_EI_SE);
q_EM_SI = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_EM_SE, q_EM_SE);
Q_EM_SI = [q_EM_SI(:,1:3).*lstarSE, q_EM_SI(:,4:6).*lstarSE./tstarSE];

t_SE = e_dep_trans/tstarSE+transfer.departureManifoldArc2.t;
q_SE = [transfer.departureManifoldArc2.x, transfer.departureManifoldArc2.y, transfer.departureManifoldArc2.z, transfer.departureManifoldArc2.xdot, transfer.departureManifoldArc2.ydot, transfer.departureManifoldArc2.zdot];
q_SE_SI = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_SE, q_SE);
Q_SE_SI = [q_SE_SI(:,1:3).*lstarSE, q_SE_SI(:,4:6).*lstarSE./tstarSE];

ode = @(t,r) ODE_C(t, r, gmS);
odeOpts = odeset('RelTol', 1E-12, 'AbsTol', 1E-12);
tspan_dep = e_dep_SoI+[0:(3600*24):transfer.departureConic.TOF, transfer.departureConic.TOF];
IC_dep = transfer.departureConic.state;
dep_sol = ode89(ode, tspan_dep, IC_dep, odeOpts);

tspan_bridge = e_bridge+[0:(3600*24):transfer.bridgeConic.TOF, transfer.bridgeConic.TOF];
IC_bridge = transfer.bridgeConic.state;
bridge_sol = ode89(ode, tspan_bridge, IC_bridge, odeOpts);

tspan_arr = e_arr_int+[0:(3600*24):transfer.arrivalConic.TOF, transfer.arrivalConic.TOF];
IC_arr = transfer.arrivalConic.state;
arr_sol = ode89(ode, tspan_arr, IC_arr, odeOpts);

t_SM = transfer.theta_arr_f+transfer.arrivalManifoldArc.t; % Because Mars is not yet realistically phased
q_SM = [transfer.arrivalManifoldArc.x, transfer.arrivalManifoldArc.y, transfer.arrivalManifoldArc.z, transfer.arrivalManifoldArc.xdot, transfer.arrivalManifoldArc.ydot, transfer.arrivalManifoldArc.zdot];
q_SM_SI = planetRotToSunEclipJ2000(muSM, transfer.initialEpoch, 'Mars', lstarSM, tstarSM, t_SM, q_SM);
Q_SM_SI = [q_SM_SI(:,1:3).*lstarSM, q_SM_SI(:,4:6).*lstarSM./tstarSM];

t_m_EM = linspace(transfer.t_0/tstarEM, e_dep_SoI/tstarEM, 10001);
q_m_EI_EM = rotToP1EclipJ2000(muEM, transfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, t_m_EM, ones(10001, 1)*[1-muEM, 0, 0, 0, 0, 0]);
q_m_EI_SE = [q_m_EI_EM(:,1:3).*lstarEM./lstarSE, q_m_EI_EM(:,4:6).*lstarEM.*tstarSE./lstarSE./tstarEM];
t_m_SE = t_m_EM.*tstarEM./tstarSE;
q_m_SE = P2EclipJ2000ToRot(muSE, transfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, t_m_SE, q_m_EI_SE);

t_E_SE = linspace(transfer.t_0/tstarSE, transfer.t_0/tstarSE+2*pi, 10001);
q_E_SI = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_E_SE, ones(10001, 1)*[1-muSE, 0, 0, 0, 0, 0]);
Q_E_SI = [q_E_SI(:,1:3).*lstarSE, q_E_SI(:,4:6).*lstarSE./tstarSE];

t_M_SM = linspace(transfer.theta_arr_f-2*pi, transfer.theta_arr_f, 10001);
q_M_SI = planetRotToSunEclipJ2000(muSM, transfer.initialEpoch, 'Mars', lstarSM, tstarSM, t_M_SM, ones(10001, 1)*[1-muSM, 0, 0, 0, 0, 0]);
Q_M_SI = [q_M_SI(:,1:3).*lstarSM, q_M_SI(:,4:6).*lstarSM./tstarSM];

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
% set(p12, 'DisplayName', "Dep. CR3BPArc")
% plot3(1-muEM+RSoIm/lstarEM*cos(linspace(0, 2*pi, 101)), RSoIm/lstarEM*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'k:', 'DisplayName', "Moon SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
% zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
% % title("Earth-Moon Rot.", 'Interpreter', 'latex')
% leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg1.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig1, 'MarsMMATCR3BP_1.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Sun-Earth Trajectory Plot
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", 10*RE/lstarSE, [1-muSE, 0, 0]);
% set(Earth, 'DisplayName', "Earth (x10)")
% p21 = plot3WithArrows(q_m_SE(:,1), q_m_SE(:,2), q_m_SE(:,3), 'k--');
% set(p21, 'DisplayName', "Moon Traj.")
% % scatter3(a1SE, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "SE $L_{1}$")
% scatter3(a2SE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "SE $L_{2}$")
% p22 = plot3WithArrows(q_EM_SE(:,1), q_EM_SE(:,2), q_EM_SE(:,3), 'r', 'NumArrows', 2, 'ArrowScale', 2);
% set(p22, 'DisplayName', "Dep. CR3BP Arc")
% p23 = plot3WithArrows(transfer.departureManifoldArc2.x(1:end-7), transfer.departureManifoldArc2.y(1:end-7), transfer.departureManifoldArc2.z(1:end-7), 'r', 'ArrowScale', 0.5);
% set(p23, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% % plot3(1-muSE+RSoIE/lstarSE*cos(linspace(0, 2*pi, 101)), RSoIE/lstarSE*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'k:', 'DisplayName', "Earth SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [SE ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [SE ndim]", 'Interpreter', 'latex')
% hz2 = zlabel("$z$ [SE ndim]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz2)
% % title("Sun-Earth Rot.", 'Interpreter', 'latex')
% leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg2.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig2, 'MarsMMATCR3BP_2.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Heliocentric MMAT Plot
% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% Sun = plot3DBody("Sun", 10*RS, [0, 0, 0]);
% set(Sun, 'DisplayName', "Sun (x10)")
% p31 = plot3WithArrows(Q_E_SI(:,1), Q_E_SI(:,2), Q_E_SI(:,3), 'g--', 'NumArrows', 2);
% set(p31, 'DisplayName', "Earth Orbit")
% p32 = plot3WithArrows(Q_M_SI(:,1), Q_M_SI(:,2), Q_M_SI(:,3), 'r--', 'NumArrows', 2, 'ArrowScale', 0.75);
% set(p32, 'DisplayName', "Mars Orbit")
% Earth = plot3DBody("Earth", 1000*RE, Q_E_SI(1,1:3));
% set(Earth, 'DisplayName', "Earth (x1000) at Dep.")
% Mars = plot3DBody("Mars", 1000*RM, Q_M_SI(end,1:3));
% set(Mars, 'DisplayName', "Mars (x1000) at Arr.")
% p33 = plot3WithArrows(Q_EM_SI(:,1), Q_EM_SI(:,2), Q_EM_SI(:,3), 'r', 'NumArrows', 2);
% set(p33, 'DisplayName', "Dep. CR3BP Arc")
% p34 = plot3WithArrows(Q_SE_SI(:,1), Q_SE_SI(:,2), Q_SE_SI(:,3), 'r');
% set(p34, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% p35 = plot3WithArrows(dep_sol.y(1,:), dep_sol.y(2,:), dep_sol.y(3,:), 'm');
% set(p35, 'DisplayName', "Dep. Conic")
% p36 = plot3WithArrows(bridge_sol.y(1,:), bridge_sol.y(2,:), bridge_sol.y(3,:), 'Color', [0.5 0 0.5]);
% set(p36, 'DisplayName', "Bridge Conic")
% p37 = plot3WithArrows(arr_sol.y(1,:), arr_sol.y(2,:), arr_sol.y(3,:), 'c', 'ArrowScale', 0.75);
% set(p37, 'DisplayName', "Arr. Conic")
% p38 = plot3WithArrows(Q_SM_SI(:,1), Q_SM_SI(:,2), Q_SM_SI(:,3), 'b', 'FlipDir', 'On', 'ArrowScale', 0.75);
% set(p38, 'DisplayName', "Arr. CR3BP Arc")
% scatter3(bridge_sol.y(1,1), bridge_sol.y(2,1), bridge_sol.y(3,1), 75, 'k', 'filled', 's', 'DisplayName', "$\Delta v_{1}="+num2str(transfer.Deltav_1)+"$ km/s")
% scatter3(bridge_sol.y(1,end), bridge_sol.y(2,end), bridge_sol.y(3,end), 75, 'k', 'filled', '^', 'DisplayName', "$\Delta v_{2}="+num2str(transfer.Deltav_2)+"$ km/s")
% axis equal
% grid on
% xlabel("$X$ [km]", 'Interpreter', 'latex')
% ylabel("$Y$ [km]", 'Interpreter', 'latex')
% hz3 = zlabel("$Z$ [km]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz3)
% % title("Sun-Centered Ecliptic J2000 $|$ TOF = "+num2str(transfer.TOF/24/3600/365.25, 3)+" yrs", 'Interpreter', 'latex')
% disp("TOF: "+transfer.TOF/24/3600/365.25+" yrs")
% leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg3.EntryContainer.NodeChildren(end-3).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig3, 'MarsMMATCR3BP_3.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Sun-Mars Trajectory Plot
% fig4 = figure("Position", [200 100 1200 750]);
% hold on
% Mars = plot3DBody("Mars", 10*RM/lstarSM, [1-muSM, 0, 0]);
% set(Mars, 'DisplayName', "Mars (x10)")
% scatter3(a1SM, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "SM $L_{1}$")
% % scatter3(a2SM, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "SM $L_{2}$")
% p41 = plot3WithArrows(transfer.arrivalOrbit.x, transfer.arrivalOrbit.y, transfer.arrivalOrbit.z, 'g--', 'NumArrows', 2, 'ArrowScale', 1.5);
% set(p41, 'DisplayName', "Arr. Orbit")
% p42 = plot3WithArrows(transfer.arrivalManifoldArc.x(1:end-6), transfer.arrivalManifoldArc.y(1:end-6), transfer.arrivalManifoldArc.z(1:end-6), 'b', 'FlipDir', 'On', 'ArrowScale', 0.5);
% set(p42, 'DisplayName', "Arr. CR3BP Arc")
% % plot3(1-muSM+RSoIM/lstarSM*cos(linspace(0, 2*pi, 101)), RSoIM/lstarSM*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'k:', 'DisplayName', "Mars SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [SM ndim]", 'Interpreter', 'latex')
% ylabel("$y$ [SM ndim]", 'Interpreter', 'latex')
% hz4 = zlabel("$z$ [SM ndim]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz4)
% % title("Sun-Mars Rot.", 'Interpreter', 'latex')
% leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig4, 'MarsMMATCR3BP_4.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Family Plots
fig5 = figure("Position", [200 100 1200 750]);
hold on
scatter(wrapTo360(theta_E_0s), TOFs/24/3600/365.25, 20, Deltavs, 'filled', 'HandleVisibility', 'off')
xlim([0 360])
grid on
xlabel("$\theta_{E,0}$ [deg]", 'Interpreter', 'latex')
ylabel("TOF [yrs]", 'Interpreter', 'latex')
% title("MMAT Family Tradespace", 'Interpreter', 'latex')
colormap(viridis)
cb5 = colorbar;
caxis([4 8])
ylabel(cb5, "$\Delta v$ [km/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
cb5.Label.Position = cb5.Label.Position+[-2 2.1 0];
set(gca, 'Color', 'w');
view(2)
hold off
% exportgraphics(fig5, 'MarsMMATCR3BP_5.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

fig6 = figure("Position", [200 100 1200 750]);
hold on
scatter(wrapTo360((6.253075709008801*180/pi):0.5:((6.253075709008801+2*pi)*180/pi)), wrapTo360((0.027114456425096738*180/pi):(0.5*tstarSE/tstarSM):((0.027114456425096738+2*pi*tstarSE/tstarSM)*180/pi)), 5, 'w', 'filled', 'HandleVisibility', 'off')
scatter(wrapTo360(theta_E_0s), wrapTo360(theta_M_0s), 20, Deltavs, 'filled', 'HandleVisibility', 'off')
axis equal
axis([0 360 0 360])
grid on
xlabel("$\theta_{E,0}$ [deg]", 'Interpreter', 'latex')
ylabel("$\theta_{M,0}$ [deg]", 'Interpreter', 'latex')
title("MMAT Family Tradespace", 'Interpreter', 'latex')
colormap(viridis)
cb6 = colorbar;
caxis([4 8])
ylabel(cb6, "$\Delta v$ [km/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
cb6.Label.Position = cb6.Label.Position+[-2 2.1 0];
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig6, 'MarsMMATCR3BP_6.png', 'BackgroundColor', 'k')

fig7 = figure("Position", [200 100 1200 750]);
hold on
scatter(TOFs/24/3600/365.25, Deltavs, 10, 'filled', 'DisplayName', "$JC_{EM}=3.03$")
yline(5.639, 'w', 'DisplayName', "Modified Hohmann Transfer") % L1
axis([2 7 4 8])
grid on
xlabel("TOF [yrs]", 'Interpreter', 'latex')
ylabel("$\Delta v$ [km/s]", 'Interpreter', 'latex')
title("$L_{1}$ Halo - $3.0001857$ $L_{1}$ Halo", 'Interpreter', 'latex')
leg7 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig7, 'MarsMMATCR3BP_7.png', 'BackgroundColor', 'k')