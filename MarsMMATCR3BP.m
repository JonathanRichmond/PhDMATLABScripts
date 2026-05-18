%%% MarsMMATCR3BP.jl
%%% Jonathan LeFevre Richmond
%%% C: 31 January 2026
%%% U: 13 May 2026

clear

%NEED TO FIX TRANSFER ANGLES
%% Import Transfer Data
transferData = load('E:/MMATData/MarsMMATCR3BP_3_13_L2Halo_flipfalse_3_000187_L1Halo_flipfalse.mat'); % Favorable-L2Lyapunov, Methodology-L2Halo
transfers = cellfun(@(n) transferData.(n), fieldnames(transferData));
n_transfers = length(transfers);
disp("Total transfers: "+n_transfers)
TOFs = zeros(n_transfers,1);
Deltav_1s = zeros(n_transfers,1);
Deltav_2s = zeros(n_transfers,1);
theta_E_deps = zeros(n_transfers,1);
theta_E_0s = zeros(n_transfers, 1);
for j = 1:n_transfers
    currentTransfer = transfers(j);
    TOFs(j) = currentTransfer.TOF;
    Deltav_1s(j) = currentTransfer.Deltav_1;
    Deltav_2s(j) = currentTransfer.Deltav_2;
    theta_E_deps(j) = wrapTo360(atan2d(currentTransfer.departureConic.state(2), currentTransfer.departureConic.state(1)));
    theta_E_0s(j) = currentTransfer.theta_dep_0;
end
Deltavs = Deltav_1s+Deltav_2s;
[~, sortDeltavs]= sort(Deltavs);
[~, sortTOFs] = sort(TOFs);
disp("Min. Deltav: "+Deltavs(sortDeltavs(1))+" km/s")
disp("Min. TOF: "+TOFs(sortTOFs(1))/24/3600/365.25+" yrs")
transferIdx = sortTOFs(100); % Favorable-sortDeltavs(85), Methodology-sortTOFs(100)
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

%% Transfer Angles
% theta_E_0s = zeros(n_transfers,1);
% theta_E_deps = zeros(n_transfers,1);
% theta_M_fs = zeros(n_transfers,1);
% transferAnglesA = zeros(n_transfers,1);
% transferAnglesB = zeros(n_transfers,1);
% parfor j = 1:n_transfers
%     currentTransfer = transfers(j);
%     t_EM_0 = currentTransfer.t_0/tstarEM+currentTransfer.departureManifoldArc1.t(1);
%     q_EM_0 = [currentTransfer.departureManifoldArc1.x(1), currentTransfer.departureManifoldArc1.y(1), currentTransfer.departureManifoldArc1.z(1), currentTransfer.departureManifoldArc1.xdot(1), currentTransfer.departureManifoldArc1.ydot(1), currentTransfer.departureManifoldArc1.zdot(1)];
%     q_EI_EM_0 = rotToP1EclipJ2000(muEM, currentTransfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, t_EM_0, q_EM_0);
%     q_EI_SE_0 = [q_EI_EM_0(:,1:3)*lstarEM/lstarSE, q_EI_EM_0(:,4:6)*lstarEM*tstarSE/lstarSE/tstarEM];
%     t_EM_SE_0 = t_EM_0*tstarEM/tstarSE;
%     q_EM_SE_0 = P2EclipJ2000ToRot(muSE, currentTransfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, t_EM_SE_0, q_EI_SE_0);
%     q_EM_SI_0 = planetRotToSunEclipJ2000(muSE, currentTransfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_EM_SE_0, q_EM_SE_0);
%     theta_E_0s(j) = wrapTo360(atan2d(q_EM_SI_0(1,2), q_EM_SI_0(1,1)));
%     theta_E_deps(j) = wrapTo360(atan2d(currentTransfer.departureConic.state(2), currentTransfer.departureConic.state(1)));
%     t_SM_0 = currentTransfer.theta_arr_f+currentTransfer.arrivalManifoldArc.t(1);
%     q_SM_0 = [currentTransfer.arrivalManifoldArc.x(1), currentTransfer.arrivalManifoldArc.y(1), currentTransfer.arrivalManifoldArc.z(1), currentTransfer.arrivalManifoldArc.xdot(1), currentTransfer.arrivalManifoldArc.ydot(1), currentTransfer.arrivalManifoldArc.zdot(1)];
%     q_SM_SI_0 = planetRotToSunEclipJ2000(muSM, currentTransfer.initialEpoch, 'Mars', lstarSM, tstarSM, t_SM_0, q_SM_0);
%     theta_M_fs(j) = wrapTo360(atan2d(q_SM_SI_0(1,2), q_SM_SI_0(1,1)));
%     transferAnglesA(j) = wrapTo360(theta_M_fs(j)-theta_E_0s(j));
%     transferAnglesB(j) = wrapTo360(theta_M_fs(j)-theta_E_deps(j));
% end

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

t_SM = transfer.theta_arr_f+transfer.arrivalManifoldArc.t; % Because Mars is not yet realistically phased
q_SM = [transfer.arrivalManifoldArc.x, transfer.arrivalManifoldArc.y, transfer.arrivalManifoldArc.z, transfer.arrivalManifoldArc.xdot, transfer.arrivalManifoldArc.ydot, transfer.arrivalManifoldArc.zdot];
q_SM_SI = planetRotToSunEclipJ2000(muSM, transfer.initialEpoch, 'Mars', lstarSM, tstarSM, t_SM, q_SM);
Q_SM_SI = [q_SM_SI(:,1:3).*lstarSM, q_SM_SI(:,4:6).*lstarSM./tstarSM];
r_SM_SI_R = (R*Q_SM_SI(:,1:3)')';

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

%% Interpolate Trajectories
odeCR3BPEM = @(t,r) ODE_CR3BP(t, r, muEM);
odeCR3BPSE = @(t,r) ODE_CR3BP(t, r, muSE);
odeCR3BPSM = @(t,r) ODE_CR3BP(t, r, muSM);

t_interp = [t_EM(1)*tstarEM:(0.5*24*3600):tspan_arr(end)-transfer.arrivalManifoldArc.TOF*tstarSM, tspan_arr(end)-transfer.arrivalManifoldArc.TOF*tstarSM];
t_dmEM_EM = [t_interp(t_interp < t_EM(end)*tstarEM)./tstarEM, t_EM(end)];
t_dmSE_SE = [t_SE(1), t_interp((t_interp > t_SE(1)*tstarSE) & (t_interp < t_SE(end)*tstarSE))./tstarSE, t_SE(end)];
t_dc_S = [tspan_dep(1), t_interp((t_interp > tspan_dep(1)) & (t_interp < tspan_dep(end))), tspan_dep(end)];
t_bc_S = [tspan_bridge(1), t_interp((t_interp > tspan_bridge(1)) & (t_interp < tspan_bridge(end))), tspan_bridge(end)];
t_ac_S = [tspan_arr(1), t_interp((t_interp > tspan_arr(1)) & (t_interp < tspan_arr(end))), tspan_arr(end)];
t_am_SM = [tspan_arr(end)/tstarSM, t_interp((t_interp > tspan_arr(end)) & (t_interp < t_interp(end)))./tstarSM, t_interp(end)/tstarSM];
t_all = [t_dmEM_EM*tstarEM, t_dmSE_SE*tstarSE, t_dc_S, t_bc_S, t_ac_S, t_am_SM*tstarSM];

IC_dmEM_EM = [transfer.departureManifoldArc1.x(1); transfer.departureManifoldArc1.y(1); transfer.departureManifoldArc1.z(1); transfer.departureManifoldArc1.xdot(1); transfer.departureManifoldArc1.ydot(1); transfer.departureManifoldArc1.zdot(1)];
sol_dmEM_EM = ode89(odeCR3BPEM, [t_EM(1) t_EM(end)], IC_dmEM_EM, odeOpts);
q_dmEM_EM = deval(sol_dmEM_EM, t_dmEM_EM)';
t_dmEM_SE = t_dmEM_EM.*tstarEM./tstarSE;
q_dmEM_EI = rotToP1EclipJ2000(muEM, transfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, t_dmEM_EM, q_dmEM_EM);
Q_dmEM_EI = [q_dmEM_EI(:,1:3).*lstarEM./lstarSE, q_dmEM_EI(:,4:6).*lstarEM.*tstarSE./lstarSE./tstarEM];
q_dmEM_SE = P2EclipJ2000ToRot(muSE, transfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, t_dmEM_SE, Q_dmEM_EI);
q_dmEM_S = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_dmEM_SE, q_dmEM_SE);
Q_dmEM_S = [q_dmEM_S(:,1:3).*lstarSE, q_dmEM_S(:,4:6).*lstarSE./tstarSE];

IC_dmSE_SE = [transfer.departureManifoldArc2.x(1); transfer.departureManifoldArc2.y(1); transfer.departureManifoldArc2.z(1); transfer.departureManifoldArc2.xdot(1); transfer.departureManifoldArc2.ydot(1); transfer.departureManifoldArc2.zdot(1)];
sol_dmSE_SE = ode89(odeCR3BPSE, [t_SE(1) t_SE(end)], IC_dmSE_SE, odeOpts);
q_dmSE_SE = deval(sol_dmSE_SE, t_dmSE_SE)';
q_dmSE_S = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, t_dmSE_SE, q_dmSE_SE);
Q_dmSE_S = [q_dmSE_S(:,1:3).*lstarSE, q_dmSE_S(:,4:6).*lstarSE./tstarSE];

Q_dc_S = deval(dep_sol, t_dc_S)';

Q_bc_S = deval(bridge_sol, t_bc_S)';

Q_ac_S = deval(arr_sol, t_ac_S)';

IC_am_SM = [transfer.arrivalManifoldArc.x(1); transfer.arrivalManifoldArc.y(1); transfer.arrivalManifoldArc.z(1); transfer.arrivalManifoldArc.xdot(1); transfer.arrivalManifoldArc.ydot(1); transfer.arrivalManifoldArc.zdot(1)];
sol_am_SM = ode89(odeCR3BPSM, [t_am_SM(end) t_am_SM(1)], IC_am_SM, odeOpts);
q_am_SM = deval(sol_am_SM, t_am_SM)';
q_am_S = planetRotToSunEclipJ2000(muSM, transfer.initialEpoch, 'Mars', lstarSM, tstarSM, t_am_SM-t_am_SM(1)+t_SM(end), q_am_SM);
Q_am_S = [q_am_S(:,1:3).*lstarSM, q_am_S(:,4:6).*lstarSM./tstarSM];

Moon_dm_EI = rotToP1EclipJ2000(muEM, transfer.initialEpoch, 'Earth', gmE, 'Moon', lstarEM, tstarEM, [t_dmEM_EM, t_dmSE_SE.*tstarSE./tstarEM], ones(length(t_dmEM_EM)+length(t_dmSE_SE), 1)*[1-muEM, zeros(1, 5)]);
MOON_dm_EI = [Moon_dm_EI(:,1:3).*lstarEM./lstarSE, Moon_dm_EI(:,4:6).*lstarEM.*tstarSE./lstarSE./tstarEM];
Moon_dm_SE = P2EclipJ2000ToRot(muSE, transfer.initialEpoch, 'Sun', gmS, 'Earth', lstarSE, tstarSE, [t_dmEM_SE, t_dmSE_SE], MOON_dm_EI);

Earth_dm_S = planetRotToSunEclipJ2000(muSE, transfer.initialEpoch, 'Earth', lstarSE, tstarSE, [t_dmEM_SE, t_dmSE_SE, t_dc_S./tstarSE, t_bc_S./tstarSE, t_ac_S./tstarSE, t_am_SM.*tstarSM./tstarSE], ones(length(t_dmEM_SE)+length(t_dmSE_SE)+length(t_dc_S)+length(t_bc_S)+length(t_ac_S)+length(t_am_SM), 1)*[1-muSE, zeros(1, 5)]);
EARTH_dm_S = [Earth_dm_S(:,1:3).*lstarSE, Earth_dm_S(:,4:6).*lstarSE./tstarSE];

t_dm_SM = ([t_dmEM_SE, t_dmSE_SE, t_dc_S./tstarSE, t_bc_S./tstarSE, t_ac_S./tstarSE, t_am_SM.*tstarSM./tstarSE]-t_dmEM_SE(1)).*tstarSE./tstarSM+transfer.theta_arr_0;
Mars_dm_S = planetRotToSunEclipJ2000(muSM, transfer.initialEpoch, 'Mars', lstarSM, tstarSM, t_dm_SM, ones(length(t_dm_SM), 1)*[1-muSM, 0, 0, 0, 0, 0]);
MARS_dm_S = [Mars_dm_S(:,1:3).*lstarSM, Mars_dm_S(:,4:6).*lstarSM./tstarSM];

%% Earth-Moon Trajectory Plot
% fig1 = figure("Position", [200 100 1200 750]);
% hold on
% % Earth = plot3DBody("Earth", RE/lstarEM, [-muEM, 0, 0].*lstarEM./lstarSE);
% % set(Earth, 'DisplayName', "Earth")
% Moon = plot3DBody("Moon", Rm/lstarEM.*lstarEM./lstarSE, [1-muEM, 0, 0].*lstarEM./lstarSE);
% set(Moon, 'DisplayName', "Moon")
% scatter3(a1EM.*lstarEM./lstarSE, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "EM $L_{1}$")
% scatter3(a2EM.*lstarEM./lstarSE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "EM $L_{2}$")
% p11 = plot3WithArrows(transfer.departureOrbit.x.*lstarEM./lstarSE, transfer.departureOrbit.y.*lstarEM./lstarSE, transfer.departureOrbit.z.*lstarEM./lstarSE, 'g--', 'NumArrows', 2, 'ArrowScale', 1.5);
% set(p11, 'DisplayName', "Dep. Orbit")
% % p12 = plot3WithArrows(transfer.departureManifoldArc1.x.*lstarEM./lstarSE, transfer.departureManifoldArc1.y.*lstarEM./lstarSE, transfer.departureManifoldArc1.z.*lstarEM./lstarSE, 'r');
% % set(p12, 'DisplayName', "Dep. CR3BP Arc")
% % plot3(1-muEM+RSoIm/lstarEM*cos(linspace(0, 2*pi, 101)), RSoIm/lstarEM*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Moon SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [AU]", 'Interpreter', 'latex')
% ylabel("$y$ [AU]", 'Interpreter', 'latex')
% zlabel("$z$ [AU]", 'Interpreter', 'latex')
% % title("Earth-Moon Rot.", 'Interpreter', 'latex')
% % leg1 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg1.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig1, 'MarsMMATCR3BP_1.png','BackgroundColor', 'w')
% % exportgraphics(fig1, 'MarsMMATCR3BP_1.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Sun-Earth Trajectory Plot
% fig2 = figure("Position", [200 100 1200 750]);
% hold on
% Earth = plot3DBody("Earth", 10*RE/lstarSE, [1-muSE, 0, 0]);
% set(Earth, 'DisplayName', "Earth")
% p21 = plot3WithArrows(q_m_SE(:,1), q_m_SE(:,2), q_m_SE(:,3), 'k--');
% set(p21, 'DisplayName', "Moon Traj.")
% scatter3(a1SE, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "SE $L_{1}$")
% % scatter3(a2SE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "SE $L_{2}$")
% p22 = plot3WithArrows(q_EM_SE(:,1), q_EM_SE(:,2), q_EM_SE(:,3), 'r', 'NumArrows', 2, 'ArrowScale', 2);
% set(p22, 'DisplayName', "Dep. CR3BP Arc")
% p23 = plot3WithArrows(transfer.departureManifoldArc2.x(1:end-9), transfer.departureManifoldArc2.y(1:end-9), transfer.departureManifoldArc2.z(1:end-9), 'r');
% set(p23, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% % plot3(1-muSE+RSoIE/lstarSE*cos(linspace(0, 2*pi, 101)), RSoIE/lstarSE*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'w:', 'DisplayName', "Earth SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [AU]", 'Interpreter', 'latex')
% ylabel("$y$ [AU]", 'Interpreter', 'latex')
% hz2 = zlabel("$z$ [AU]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz2)
% % title("Sun-Earth Rot.", 'Interpreter', 'latex')
% % leg2 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg2.EntryContainer.NodeChildren(end).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig2, 'MarsMMATCR3BP_2.png', 'BackgroundColor', 'w')
% % exportgraphics(fig2, 'MarsMMATCR3BP_2.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Heliocentric MMAT Plot
disp("TOF: "+transfer.TOF/24/3600/365.25+" yrs")

% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% Sun = plot3DBody("Sun", 10*RS./lstarSE, [0, 0, 0]);
% set(Sun, 'DisplayName', "Sun")
% p31 = plot3WithArrows(Q_E_SI(:,1)./lstarSE, Q_E_SI(:,2)./lstarSE, Q_E_SI(:,3)./lstarSE, 'g--', 'NumArrows', 2);
% set(p31, 'DisplayName', "Earth Orbit")
% p32 = plot3WithArrows(Q_M_SI(:,1)./lstarSE, Q_M_SI(:,2)./lstarSE, Q_M_SI(:,3)./lstarSE, 'r--', 'NumArrows', 2, 'ArrowScale', 0.75);
% set(p32, 'DisplayName', "Mars Orbit")
% Earth = plot3DBody("Earth", 1000*RE./lstarSE, Q_E_SI(1,1:3)./lstarSE);
% set(Earth, 'DisplayName', "Earth at Dep.")
% Mars = plot3DBody("Mars", 1000*RM./lstarSE, Q_M_SI(end,1:3)./lstarSE);
% set(Mars, 'DisplayName', "Mars at Arr.")
% p33 = plot3WithArrows(Q_EM_SI(:,1)./lstarSE, Q_EM_SI(:,2)./lstarSE, Q_EM_SI(:,3)./lstarSE, 'r', 'NumArrows', 1, 'ArrowScale', 10);
% set(p33, 'DisplayName', "Dep. CR3BP Arc")
% p34 = plot3WithArrows(Q_SE_SI(:,1)./lstarSE, Q_SE_SI(:,2)./lstarSE, Q_SE_SI(:,3)./lstarSE, 'r', 'NumArrows', 2, 'ArrowScale', 1.5);
% set(p34, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% p35 = plot3WithArrows(dep_sol.y(1,:)./lstarSE, dep_sol.y(2,:)./lstarSE, dep_sol.y(3,:)./lstarSE, 'm', 'NumArrows', 1, 'ArrowScale', 4);
% set(p35, 'DisplayName', "Dep. Conic")
% p36 = plot3WithArrows(bridge_sol.y(1,:)./lstarSE, bridge_sol.y(2,:)./lstarSE, bridge_sol.y(3,:)./lstarSE, 'Color', [0.5 0 0.5]);
% set(p36, 'DisplayName', "Bridge Conic")
% p37 = plot3WithArrows(arr_sol.y(1,:)./lstarSE, arr_sol.y(2,:)./lstarSE, arr_sol.y(3,:)./lstarSE, 'c', 'NumArrows', 1, 'ArrowScale', 15);
% set(p37, 'DisplayName', "Arr. Conic")
% p38 = plot3WithArrows(Q_SM_SI(:,1)./lstarSE, Q_SM_SI(:,2)./lstarSE, Q_SM_SI(:,3)./lstarSE, 'b', 'FlipDir', 'On', 'ArrowScale', 0.75);
% set(p38, 'DisplayName', "Arr. CR3BP Arc")
% scatter3(bridge_sol.y(1,1)./lstarSE, bridge_sol.y(2,1)./lstarSE, bridge_sol.y(3,1)./lstarSE, 75, 'k', 'filled', 's', 'DisplayName', "$\Delta v_{1}="+num2str(transfer.Deltav_1)+"$ km/s")
% scatter3(bridge_sol.y(1,end)./lstarSE, bridge_sol.y(2,end)./lstarSE, bridge_sol.y(3,end)./lstarSE, 75, 'k', 'filled', '^', 'DisplayName', "$\Delta v_{2}="+num2str(transfer.Deltav_2)+"$ km/s")
% axis equal
% grid on
% xlabel("$X$ [AU]", 'Interpreter', 'latex')
% ylabel("$Y$ [AU]", 'Interpreter', 'latex')
% hz3 = zlabel("$Z$ [AU]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz3)
% % title("Sun-Centered Ecliptic J2000 $|$ TOF = "+num2str(transfer.TOF/24/3600/365.25, 3)+" yrs", 'Interpreter', 'latex')
% % leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % drawnow;
% % set(leg3.EntryContainer.NodeChildren(end-3).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(3)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig3, 'MarsMMATCR3BP_3.png', 'BackgroundColor', 'w')
% % exportgraphics(fig3, 'MarsMMATCR3BP_3.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

% fig3 = figure("Position", [200 100 1200 750]);
% hold on
% Sun = plot3DBody("Sun", 10*RS, [0, 0, 0]);
% set(Sun, 'DisplayName', "Sun (x10)")
% p31 = plot3WithArrows(Q_E_SI(:,1), Q_E_SI(:,2), Q_E_SI(:,3), 'g--', 'NumArrows', 2);
% set(p31, 'DisplayName', "Earth Orbit")
% p32 = plot3WithArrows(Q_M_SI(:,1), Q_M_SI(:,2), Q_M_SI(:,3), 'r--', 'NumArrows', 2, 'ArrowScale', 0.75);
% set(p32, 'DisplayName', "Mars Orbit")
% Earth = plot3DBody("Earth", 1000*RE, [lstarSE, 0, 0]);
% set(Earth, 'DisplayName', "Earth (x1000) at Dep. Conic")
% Mars = plot3DBody("Mars", 1000*RM, (R*Q_M_SI(end,1:3)')');
% set(Mars, 'DisplayName', "Mars (x1000) at Arr.")
% p33 = plot3WithArrows(r_EM_SI_R(:,1), r_EM_SI_R(:,2), r_EM_SI_R(:,3), 'r', 'NumArrows', 2);
% set(p33, 'DisplayName', "Dep. CR3BP Arc")
% p34 = plot3WithArrows(r_SE_SI_R(:,1), r_SE_SI_R(:,2), r_SE_SI_R(:,3), 'r');
% set(p34, 'DisplayName', "Dep. CR3BP Arc", 'HandleVisibility', 'off')
% p35 = plot3WithArrows(r_dep_R(1,:), r_dep_R(2,:), r_dep_R(3,:), 'm');
% set(p35, 'DisplayName', "Dep. Conic Rot.")
% p36 = plot3WithArrows(r_bridge_R(1,:), r_bridge_R(2,:), r_bridge_R(3,:), 'Color', [0.5 0 0.5]);
% set(p36, 'DisplayName', "Bridge Conic")
% p37 = plot3WithArrows(r_arr_R(1,:), r_arr_R(2,:), r_arr_R(3,:), 'c', 'ArrowScale', 0.75);
% set(p37, 'DisplayName', "Arr. Conic")
% p38 = plot3WithArrows(r_SM_SI_R(:,1), r_SM_SI_R(:,2), r_SM_SI_R(:,3), 'b', 'FlipDir', 'On', 'ArrowScale', 0.75);
% set(p38, 'DisplayName', "Arr. CR3BP Arc")
% scatter3(r_bridge_R(1,1), r_bridge_R(2,1), r_bridge_R(3,1), 75, 'k', 'filled', 's', 'DisplayName', "$\Delta v_{1}="+num2str(transfer.Deltav_1)+"$ km/s")
% scatter3(r_bridge_R(1,end), r_bridge_R(2,end), r_bridge_R(3,end), 75, 'k', 'filled', '^', 'DisplayName', "$\Delta v_{2}="+num2str(transfer.Deltav_2)+"$ km/s")
% axis equal
% grid on
% xlabel("$X$ [km]", 'Interpreter', 'latex')
% ylabel("$Y$ [km]", 'Interpreter', 'latex')
% title("SCI Ecliptic J2000 $|$ TOF = "+num2str(transfer.TOF/24/3600/365.25, 3)+" yrs", 'Interpreter', 'latex')
% % title("$\phi_{H}$ = "+num2str(wrapTo180(transferAnglesB(transferIdx)), 4)+" deg", 'Interpreter', 'latex')
% leg3 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% drawnow;
% set(leg3.EntryContainer.NodeChildren(end-3).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]))
% set(gca, 'Color', 'w');
% view(2)
% hold off
% % ax = gca;
% % ax.SortMethod = 'childorder';
% % exportgraphics(fig3, 'MarsMMATCR3BP_3.png', 'BackgroundColor', 'w')

%% Sun-Mars Trajectory Plot
% fig4 = figure("Position", [200 100 1200 750]);
% hold on
% Mars = plot3DBody("Mars", 10*RM/lstarSM.*lstarSM./lstarSE, [1-muSM, 0, 0].*lstarSM./lstarSE);
% set(Mars, 'DisplayName', "Mars")
% scatter3(a1SM.*lstarSM./lstarSE, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "SM $L_{1}$")
% % scatter3(a2SM.*lstarSM./lstarSE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "SM $L_{2}$")
% p41 = plot3WithArrows(transfer.arrivalOrbit.x.*lstarSM./lstarSE, transfer.arrivalOrbit.y.*lstarSM./lstarSE, transfer.arrivalOrbit.z.*lstarSM./lstarSE, 'g--', 'NumArrows', 2, 'ArrowScale', 1.5);
% set(p41, 'DisplayName', "Arr. Orbit")
% p42 = plot3WithArrows(q_am_SM(101:end,1).*lstarSM./lstarSE, q_am_SM(101:end,2).*lstarSM./lstarSE, q_am_SM(101:end,3).*lstarSM./lstarSE, 'b', 'ArrowScale', 0.5);
% set(p42, 'DisplayName', "Arr. CR3BP Arc")
% % plot3(1-muSM+RSoIM/lstarSM*cos(linspace(0, 2*pi, 101)), RSoIM/lstarSM*sin(linspace(0, 2*pi, 101)), zeros(1, 101), 'k:', 'DisplayName', "Mars SoI Radius")
% axis equal
% grid on
% xlabel("$x$ [AU]", 'Interpreter', 'latex')
% ylabel("$y$ [AU]", 'Interpreter', 'latex')
% hz4 = zlabel("$z$ [AU]", 'Interpreter', 'latex');
% shiftZLabel(gca, hz4, 0.1)
% % title("Sun-Mars Rot.", 'Interpreter', 'latex')
% % leg4 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% set(gca, 'Color', 'w');
% view(3)
% hold off
% ax = gca;
% ax.SortMethod = 'childorder';
% % exportgraphics(fig4, 'MarsMMATCR3BP_4.png', 'BackgroundColor', 'w')
% % exportgraphics(fig4, 'MarsMMATCR3BP_4.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Family Plots
years = 0; % 2030
DeltavHohmann = 5.694; % km/s (Mars)

% fig5 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(wrapTo2Pi(theta_E_0s).*180./pi, TOFs/24/3600/365.25, 20, Deltavs, 'filled', 'HandleVisibility', 'off')
% xlim([0 360])
% grid on
% xlabel("$\theta_{E,0}$ [deg]", 'Interpreter', 'latex')
% ylabel("TOF [yrs]", 'Interpreter', 'latex')
% % title("MMAT Family Tradespace", 'Interpreter', 'latex')
% colormap(viridis)
% cb5 = colorbar;
% caxis([4 8])
% ylabel(cb5, "$\Delta v$ [km/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb5.Label.Position = cb5.Label.Position+[-2 2.1 0];
% set(gca, 'Color', 'w');
% view(2)
% hold off
% % exportgraphics(fig5, 'MarsMMATCR3BP_5.png', 'BackgroundColor', 'w')

% fig6 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(wrapTo360(((6.253075709008801+years*2*pi)*180/pi):0.5:((6.253075709008801+(years+1)*2*pi)*180/pi)), wrapTo360(((0.027114456425096738+years*2*pi*tstarSE/tstarSM)*180/pi):(0.5*tstarSE/tstarSM):((0.027114456425096738+(years+1)*2*pi*tstarSE/tstarSM)*180/pi)), 5, 'w', 'filled', 'HandleVisibility', 'off')
% scatter(wrapTo360(theta_E_0s), wrapTo360(theta_M_0s), 20, Deltavs, 'filled', 'HandleVisibility', 'off')
% axis equal
% axis([0 360 0 360])
% grid on
% xlabel("$\theta_{E,0}$ [deg]", 'Interpreter', 'latex')
% ylabel("$\theta_{M,0}$ [deg]", 'Interpreter', 'latex')
% % title("MMAT Family Tradespace", 'Interpreter', 'latex')
% colormap(viridis)
% cb6 = colorbar;
% caxis([4 8])
% ylabel(cb6, "$\Delta v$ [km/s]", 'Interpreter', 'latex', 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% cb6.Label.Position = cb6.Label.Position+[-2 2.1 0];
% set(gca, 'Color', 'w');
% view(2)
% hold off
% % exportgraphics(fig6, 'MarsMMATCR3BP_6.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Comparison Metrics
mask = Deltavs < DeltavHohmann;
DeltavsComp = Deltavs(mask);
TOFsComp = TOFs(mask);
[sortTOFsComp, sortTOFsCompIdx] = sort(TOFsComp);
sortDeltavsComp = DeltavsComp(sortTOFsCompIdx);

dTOF = diff(sortTOFsComp);
breaks = find(dTOF > 0.025*3600*24*365.25)';
segments = [1, breaks+1; breaks, length(sortTOFsComp)];
[~, keepIndex] = max(segments(2,:)-segments(1,:));
sortTOFsKeep = sortTOFsComp(segments(1,keepIndex):segments(2,keepIndex));
sortTOFsKeepIdx = sortTOFsCompIdx(segments(1,keepIndex):segments(2,keepIndex));
sortDeltavsKeep = sortDeltavsComp(segments(1,keepIndex):segments(2,keepIndex));

sortParetoMask = false(size(sortTOFsKeep));
bestDeltav = inf;
for j = 1:length(sortTOFsKeep)
    if sortDeltavsKeep(j) < bestDeltav
        sortParetoMask(j) = true;
        bestDeltav = sortDeltavsKeep(j);
    end
end
paretoMask = sortParetoMask;
TOFsPareto = sortTOFsKeep(paretoMask);
DeltavsPareto = sortDeltavsKeep(paretoMask);
[sortTOFsPareto, sortTOFsParetoIdx] = sort(TOFsPareto/24/3600/365.25);
sortDeltavsPareto = DeltavsPareto(sortTOFsParetoIdx);
for j = 2:length(sortTOFsPareto)
    if sortDeltavsPareto(j) >= sortDeltavsPareto(j-1)
        sortDeltavsPareto(j) = sortDeltavsPareto(j-1)-1E-6;
    end
end

pLeft = polyfit(sortTOFsPareto(1:5), sortDeltavsPareto(1:5), 1);
slope0 = min(pLeft(1), -1.5);
TOFHohmann = sortTOFsPareto(2)+(DeltavHohmann-sortDeltavsPareto(2))/slope0;
pRight = polyfit(sortTOFsPareto(end-4:end), sortDeltavsPareto(end-4:end), 1);
slope1 = min(pRight(1), -1);
TOFRight = sortTOFsPareto(end-1)+(sortDeltavsPareto(end)-sortDeltavsPareto(end-1))/slope1;
pp = pchip([TOFHohmann; sortTOFsPareto(2:end-1); TOFRight], [DeltavHohmann; sortDeltavsPareto(2:end)]);
TOFSmooth = linspace(TOFHohmann, TOFRight, 500);
DeltavFixInterp = ppval(pp, TOFSmooth);
DeltavSmooth = sgolayfilt(DeltavFixInterp, 2, 201);
slopes = diff(DeltavSmooth)./diff(TOFSmooth);
startIndex = find((slopes < -1) & (DeltavSmooth(1:end-1) < DeltavHohmann), 1, 'first');
if startIndex == 1
    startIndex = startIndex+1;
end
cutIndex = find((slopes < -0.5) & (DeltavSmooth(1:end-1) > sortDeltavsPareto(end)), 1, 'last');
if cutIndex == length(TOFSmooth)
    cutIndex = cutIndex-1;
end
DeltavMono = DeltavSmooth(startIndex:cutIndex);
DeltavFine = [DeltavHohmann, DeltavMono, sortDeltavsPareto(end)];
TOFMono = TOFSmooth(startIndex:cutIndex);
pLeftNew = polyfit(TOFMono(1:2), DeltavMono(1:2), 1);
slope0New = pLeftNew(1);
TOFHohmannNew = TOFMono(1)+(DeltavHohmann-DeltavMono(1))/slope0New;
pRightNew = polyfit(TOFMono(end-1:end), DeltavMono(end-1:end), 1);
slope1New = pRightNew(1);
TOFRightNew = TOFMono(end)+(sortDeltavsPareto(end)-DeltavMono(end))/slope1New;
TOFFine = [TOFHohmannNew, TOFMono, TOFRightNew];

refTOF = 8;
HV = trapz([TOFFine, refTOF], DeltavHohmann*ones(1,length(DeltavFine)+1)-[DeltavFine, DeltavFine(end)]);
disp("Hypervolume: "+HV+" yrs*km/s")

plotShift = 0.025;
DeltavRef = DeltavHohmann*ones(length(TOFFine)+1,1);
xFill = [[TOFFine'; refTOF+plotShift]; flipud([TOFFine'; refTOF+plotShift])];
yFill = [DeltavRef; flipud([DeltavFine'; DeltavFine(end)])];

pointsMask = false(length(TOFsComp),1);
for j = 1:length(TOFsComp)
    if TOFsComp(j)/24/3600/365.25 >= interp1(DeltavFine, TOFFine, DeltavsComp(j), 'pchip')
        pointsMask(j) = true;
    end
end
TOFsCons = TOFsComp(pointsMask);
DeltavsCons = DeltavsComp(pointsMask);

TOFs_n = (TOFsCons-min(TOFsCons))/(max(TOFsCons)-min(TOFsCons));
Deltavs_n = (DeltavsCons-min(DeltavsCons))/(max(DeltavsCons)-min(DeltavsCons));
TOFsPareto_n = (TOFFine-min(TOFFine))/(max(TOFFine)-min(TOFFine));
DeltavsPareto_n = (DeltavFine-min(DeltavFine))/(max(DeltavFine)-min(DeltavFine));
pareto = [TOFsPareto_n', DeltavsPareto_n'];
points = [TOFs_n, Deltavs_n];
distPareto = zeros(size(points, 1),1);
for j = 1:size(points, 1)
    d = sqrt(sum((pareto-points(j,:)).^2, 2));
    distPareto(j) = min(d);
end
DM = mean(distPareto);
disp("Mean distance: "+DM)

data = [TOFsCons/24/3600/365.25, DeltavsCons];
xgrid = linspace(min(TOFFine), max(TOFsCons)/24/3600/365.25, 200);
ygrid = linspace(min(DeltavsCons), max(DeltavsCons), 200);
[Xg, Yg] = meshgrid(xgrid, ygrid);
gridPoints = [Xg(:), Yg(:)];
[f, ~] = ksdensity(data, gridPoints);
F = reshape(f, size(Xg));
mask = Xg < interp1(DeltavFine, TOFFine, Yg, 'pchip')-plotShift;
F(mask) = NaN;

% fig7 = figure("Position", [200 100 1200 750]);
% hold on
% scatter(TOFs/24/3600/365.25, Deltavs, 10, 'filled', 'MarkerFaceAlpha', 0.25, 'DisplayName', "Transfers")
% % scatter(TOFs/24/3600/365.25, Deltavs, 10, angleColor(deg2rad(transferAnglesA)), 'filled', 'HandleVisibility', 'off')
% % scatter(TOFs/24/3600/365.25, Deltavs, 10, angleColor(deg2rad(transferAnglesB)), 'filled', 'HandleVisibility', 'off')
% % scatter(TOFs/24/3600/365.25, Deltavs, 10, angleColor(deg2rad(transferAnglesA-transferAnglesB)), 'filled', 'HandleVisibility', 'off')
% scatter(TOFs(transferIdx)/24/3600/365.25, Deltavs(transferIdx), 30, 'k', 'filled', 'HandleVisibility', 'off')
% yline(DeltavHohmann, 'k--', 'LineWidth', 2, 'DisplayName', "Hohmann Baseline")
% plot(TOFFine-plotShift, DeltavFine, 'Color', [0.75 0 0], 'DisplayName', "Pareto Front")
% fill(xFill-plotShift, yFill, [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', "Hypervolume") % 0.1
% contour(Xg, Yg, F, 5, 'LineColor', 0.5*[0 0 0], 'LineWidth', 1, 'DisplayName', 'Density Contours')
% axis([2 8 4 9])
% grid on
% xlabel("TOF [yrs]", 'Interpreter', 'latex')
% ylabel("$\Delta v$ [km/s]", 'Interpreter', 'latex')
% % title("Earth-Moon $3.13$ $L_{2}$ S Halo - Sun-Mars $3.000187$ $L_{1}$ S Halo", 'Interpreter', 'latex')
% % leg7 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% % phasemap;
% % pb7 = phasebar('deg', 'Location', 'northeast', 'Size', 0.275, 'Title', "Test");
% % tpb7 = title(pb7, '$\phi$ [deg]', 'Interpreter', 'latex');
% % tpb7 = title(pb7, '$\phi_{H}$ [deg]', 'Interpreter', 'latex');
% % tpb7 = title(pb7, '$\phi_{M}$ [deg]', 'Interpreter', 'latex');
% % tpb7.Position(2) = tpb7.Position(2)-40;
% set(gca, 'Color', 'w');
% view(2)
% hold off
% % exportgraphics(fig7, 'MarsMMATCR3BP_7.png', 'BackgroundColor', 'w')
% % exportgraphics(fig7, 'MarsMMATCR3BP_8.pdf', 'BackgroundColor', 'w', 'ContentType', 'vector')

%% Family Comparison
% % [142:167, 169:194, 196:221]
% sheetData = readcell('Transfer Metrics.xlsx', 'Sheet', 'Mars Orientation', 'Range', 'A196:L221');
% JC_Lyap = cell2mat(sheetData(6:11,3));
% JC_halo = cell2mat(sheetData(1:4,3));
% JC_vert = cell2mat(sheetData(13:18,3));
% JC_axial = cell2mat(sheetData(20:21,3));
% JC_butt = cell2mat(sheetData(23:26,3));
% HV_Lyap = cell2mat(sheetData(6:11,11));
% HV_halo = cell2mat(sheetData(1:4,11));
% HV_vert = cell2mat(sheetData(13:18,11));
% HV_axial = cell2mat(sheetData(20:21,11));
% HV_butt = cell2mat(sheetData(23:26,11));
% MD_Lyap = cell2mat(sheetData(6:11,12));
% MD_halo = cell2mat(sheetData(1:4,12));
% MD_vert = cell2mat(sheetData(13:18,12));
% MD_axial = cell2mat(sheetData(20:21,12));
% MD_butt = cell2mat(sheetData(23:26,12));
% 
% fig9 = figure("Position", [200 100 1200 750]);
% tiledlayout(2,1)
% ax1 = nexttile;
% hold on
% plot(JC_Lyap, HV_Lyap, 'r*-', 'DisplayName', "$L_{2}$ Lyapunov")
% plot(JC_halo, HV_halo, 'g*-', 'DisplayName', "$L_{2}$ N Halo")
% plot(JC_vert, HV_vert, 'b*-', 'DisplayName', "$L_{2}$ Vertical")
% plot(JC_axial, HV_axial, 'm*-', 'DisplayName', "$L_{2}$ SW Axial")
% plot(JC_butt, HV_butt, 'c*-', 'DisplayName', "$L_{2}$ S Butterfly")
% xticks([2.97, 3.0, 3.03, 3.07, 3.1, 3.13])
% ylim([2.5 5.5])
% grid on
% ylabel("Hypervolume", 'Interpreter', 'latex')
% leg9 = legend('Location', 'bestoutside', 'Interpreter', 'latex');
% hold off
% ax2 = nexttile;
% hold on
% plot(JC_Lyap, MD_Lyap, 'r*-', 'DisplayName', "$L_{2}$ Lyapunov")
% plot(JC_halo, MD_halo, 'g*-', 'DisplayName', "$L_{2}$ S Halo")
% plot(JC_vert, MD_vert, 'b*-', 'DisplayName', "$L_{2}$ Vertical")
% plot(JC_axial, MD_axial, 'm*-', 'DisplayName', "$L_{2}$ SW Axial")
% plot(JC_butt, MD_butt, 'c*-', 'DisplayName', "$L_{2}$ S Butterfly")
% xticks([2.97, 3.0, 3.03, 3.07, 3.1, 3.13])
% ylim([0.12 0.22])
% grid on
% ylabel("Mean Distance", 'Interpreter', 'latex')
% hold off
% xlabel("Earth-Moon JC", 'Interpreter', 'latex')
% % title("Earth-Moon Departure Orbit - Sun-Mars $3.000187$ $L_{2}$ Southern Halo", 'Interpreter', 'latex')
% set(gca, 'Color', 'w');
% % exportgraphics(fig9, 'MarsMMATCR3BP_9.png', 'BackgroundColor', 'w')

%% Animation
% fig10 = figure('Position', [100 100 1500 750]);
% 
% axS = axes('Position', [0.08 0.04 0.52 0.91]);
% insetPos = [0.67 0.4 0.23 0.55];
% axEM = axes('Position', insetPos);
% axSE = axes('Position', insetPos);
% axSM = axes('Position', insetPos);
% allAxes = [axS, axEM, axSE, axSM];
% for a = allAxes
%     axis(a, 'equal');
%     grid(a, 'on');
%     hold(a, 'on');
% end
% 
% xlabel(axS, "$X$ [AU]", 'Interpreter', 'latex')
% ylabel(axS, "$Y$ [AU]", 'Interpreter', 'latex')
% zlabel(axS, "$Z$ [AU]", 'Interpreter', 'latex')
% title(axS, "Sun-Centered Eclip. J2000", 'Interpreter', 'latex')
% view(axS, 3)
% xlabel(axEM, "$x$ [AU]", 'Interpreter', 'latex')
% ylabel(axEM, "$y$ [AU]", 'Interpreter', 'latex')
% zlabel(axEM, "$z$ [AU]", 'Interpreter', 'latex')
% title(axEM, "Earth-Moon Rot.", 'Interpreter', 'latex')
% view(axEM, 3)
% xlabel(axSE, "$x$ [AU]", 'Interpreter', 'latex')
% ylabel(axSE, "$y$ [AU]", 'Interpreter', 'latex')
% hzSE = zlabel(axSE, "$z$ [AU]", 'Interpreter', 'latex');
% title(axSE, "Sun-Earth Rot.", 'Interpreter', 'latex')
% view(axSE, 3)
% xlabel(axSM, "$x$ [AU]", 'Interpreter', 'latex')
% ylabel(axSM, "$y$ [AU]", 'Interpreter', 'latex')
% zlabel(axSM, "$z$ [AU]", 'Interpreter', 'latex')
% title(axSM, "Sun-Mars Rot.", 'Interpreter', 'latex')
% view(axSM, 3)
% 
% [Sun_S, ~, ~, ~] = plot3DBodyVid(axS, "Sun", 10*RS/lstarSE, [0, 0, 0]);
% set(Sun_S, 'DisplayName', "Sun")
% theta = 0:0.01:2*pi;
% traj_Earth_S = plot3(axS, EARTH_dm_S(1:740,1)./lstarSE, EARTH_dm_S(1:740,2)./lstarSE, EARTH_dm_S(1:740,3)./lstarSE, 'g:', 'DisplayName', "Earth Orbit");
% traj_Mars_S = plot3(axS, MARS_dm_S(1:1390,1)./lstarSE, MARS_dm_S(1:1390,2)./lstarSE, MARS_dm_S(1:1390,3)./lstarSE, 'r:', 'DisplayName', "Mars Orbit");
% traj_dmEM_S = plot3(axS, Q_dmEM_S(:,1)./lstarSE, Q_dmEM_S(:,2)./lstarSE, Q_dmEM_S(:,3)./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% traj_dmSE_S = plot3(axS, Q_dmSE_S(:,1)./lstarSE, Q_dmSE_S(:,2)./lstarSE, Q_dmSE_S(:,3)./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% traj_dc_S = plot3(axS, Q_dc_S(:,1)./lstarSE, Q_dc_S(:,2)./lstarSE, Q_dc_S(:,3)./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% traj_bc_S = plot3(axS, Q_bc_S(:,1)./lstarSE, Q_bc_S(:,2)./lstarSE, Q_bc_S(:,3)./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% traj_ac_S = plot3(axS, Q_ac_S(:,1)./lstarSE, Q_ac_S(:,2)./lstarSE, Q_ac_S(:,3)./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% traj_am_S = plot3(axS, Q_am_S(:,1)./lstarSE, Q_am_S(:,2)./lstarSE, Q_am_S(:,3)./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% Deltav1 = scatter3(axS, NaN, NaN, NaN, 20, 'k', 'filled', '^', 'DisplayName', "Maneuver");
% Deltav2 = scatter3(axS, NaN, NaN, NaN, 20, 'k', 'filled', '^', 'HandleVisibility', 'off');
% 
% % [Earth_EM, ~, ~, ~] = plot3DBodyVid(axEM, "Earth", RE/lstarSE, [-muEM, 0, 0].*lstarEM./lstarSE);
% % set(Earth_EM, 'DisplayName', "Earth")
% [Moon_EM, ~, ~, ~] = plot3DBodyVid(axEM, "Moon", Rm/lstarSE, [1-muEM, 0, 0].*lstarEM./lstarSE);
% set(Moon_EM, 'DisplayName', "Moon")
% L1_EM = scatter3(axEM, a1EM*lstarEM/lstarSE, 0, 0, 20, 'r', 'filled', 'd', 'HandleVisibility', 'off');
% L2_EM = scatter3(axEM, a2EM*lstarEM/lstarSE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'HandleVisibility', 'off');
% orbit_EM = plot3(axEM, transfer.departureOrbit.x.*lstarEM./lstarSE, transfer.departureOrbit.y.*lstarEM./lstarSE, transfer.departureOrbit.z.*lstarEM./lstarSE, ':', 'Color', [0.5 0.5 0], 'DisplayName', "Dep. Orbit");
% traj_dmEM_EM = plot3(axEM, q_dmEM_EM(:,1).*lstarEM./lstarSE, q_dmEM_EM(:,2).*lstarEM./lstarSE, q_dmEM_EM(:,3).*lstarEM./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% 
% [Earth_SE, ~, ~, ~] = plot3DBodyVid(axSE, "Earth", 10*RE/lstarSE, [1-muSE, 0, 0]);
% set(Earth_SE, 'HandleVisibility', 'off')
% % L1_SE = scatter3(axSE, a1SE, 0, 0, 20, 'r', 'filled', 'd', 'HandleVisibility', 'off');
% L2_SE = scatter3(axSE, a2SE, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'HandleVisibility', 'off');
% traj_Moon_SE = plot3(axSE, Moon_dm_SE(:,1), Moon_dm_SE(:,2), Moon_dm_SE(:,3), 'k:', 'DisplayName', "Moon Orbit");
% traj_dmEM_SE = plot3(axSE, q_dmEM_SE(:,1), q_dmEM_SE(:,2), q_dmEM_SE(:,3), 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% traj_dmSE_SE = plot3(axSE, q_dmSE_SE(:,1), q_dmSE_SE(:,2), q_dmSE_SE(:,3), 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% 
% [Mars_SM, ~, ~, ~] = plot3DBodyVid(axSM, "Mars", 10*RM/lstarSE, [1-muSM, 0, 0].*lstarSM./lstarSE);
% set(Mars_SM, 'DisplayName', "Mars")
% L1_SM = scatter3(axSM, a1SM*lstarSM/lstarSE, 0, 0, 20, 'r', 'filled', 'd', 'HandleVisibility', 'off');
% orbit_SM = plot3(axSM, transfer.arrivalOrbit.x.*lstarSM./lstarSE, transfer.arrivalOrbit.y.*lstarSM./lstarSE, transfer.arrivalOrbit.z.*lstarSM./lstarSE, ':', 'Color', [0 0.5 0.5], 'DisplayName', "Arr. Orbit");
% traj_am_SM = plot3(axSM, q_am_SM(:,1).*lstarSM./lstarSE, q_am_SM(:,2).*lstarSM./lstarSE, q_am_SM(:,3).*lstarSM./lstarSE, 'Color', [0 0 0 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
% 
% hist_dmEM_S = plot3(axS, NaN, NaN, NaN, 'r', 'DisplayName', "Dep. CR3BP Arc");
% hist_dmSE_S = plot3(axS, NaN, NaN, NaN, 'r', 'HandleVisibility', 'off');
% hist_dc_S = plot3(axS, NaN, NaN, NaN, 'm', 'DisplayName', "Dep. Conic");
% hist_bc_S = plot3(axS, NaN, NaN, NaN, 'Color', [0.5 0 0.5], 'DisplayName', "Bridge Conic");
% hist_ac_S = plot3(axS, NaN, NaN, NaN, 'c', 'DisplayName', "Arr. Conic");
% hist_am_S = plot3(axS, NaN, NaN, NaN, 'b', 'DisplayName', "Arr. CR3BP Arc");
% 
% hist_dmEM_EM = plot3(axEM, NaN, NaN, NaN, 'r', 'HandleVisibility', 'off');
% 
% hist_dmEM_SE = plot3(axSE, NaN, NaN, NaN, 'r', 'HandleVisibility', 'off');
% hist_dmSE_SE = plot3(axSE, NaN, NaN, NaN, 'r', 'HandleVisibility', 'off');
% 
% hist_amSM_SM = plot3(axSM, NaN, NaN, NaN, 'b' , 'HandleVisibility', 'off');
% 
% [mark_Earth_S, sx_E, sy_E, sz_E] = plot3DBodyVid(axS, "Earth", 1000*RE/lstarSE, [0, 0, 0]);
% [mark_Mars_S, sx_M, sy_M, sz_M] = plot3DBodyVid(axS, "Mars", 1000*RM/lstarSE, [0, 0, 0]);
% mark_dm_S = scatter3(axS, NaN, NaN, NaN, 20, 'r', 'filled', 'HandleVisibility', 'off');
% mark_dc_S = scatter3(axS, NaN, NaN, NaN, 20, 'm', 'filled', 'HandleVisibility', 'off');
% mark_bc_S = scatter3(axS, NaN, NaN, NaN, 20, [0.5 0 0.5], 'filled', 'HandleVisibility', 'off');
% mark_ac_S = scatter3(axS, NaN, NaN, NaN, 20, 'c', 'filled', 'HandleVisibility', 'off');
% mark_am_S = scatter3(axS, NaN, NaN, NaN, 20, 'b', 'filled', 'HandleVisibility', 'off');
% 
% mark_EM = scatter3(axEM, NaN, NaN, NaN, 20, 'r', 'filled', 'HandleVisibility', 'off');
% 
% [mark_Moon_SE, sx_m, sy_m, sz_m] = plot3DBodyVid(axSE, "Moon", 10*Rm/lstarSE, [0, 0, 0]);
% mark_SE = scatter3(axSE, NaN, NaN, NaN, 20, 'r', 'filled', 'HandleVisibility', 'off');
% 
% mark_SM = scatter3(axSM, NaN, NaN, NaN, 20, 'b', 'filled', 'HandleVisibility', 'off');
% 
% pad = 0.05;
% allx_S = [-10*RS/lstarSE; 10*RS/lstarSE; EARTH_dm_S(:,1)./lstarSE; MARS_dm_S(:,1)./lstarSE; Q_dmEM_S(:,1)./lstarSE; Q_dmSE_S(:,1)./lstarSE; Q_dc_S(:,1)./lstarSE; Q_bc_S(:,1)./lstarSE; Q_ac_S(:,1)./lstarSE; Q_am_S(:,1)./lstarSE];
% ally_S = [-10*RS/lstarSE; 10*RS/lstarSE; EARTH_dm_S(:,2)./lstarSE; MARS_dm_S(:,2)./lstarSE; Q_dmEM_S(:,2)./lstarSE; Q_dmSE_S(:,2)./lstarSE; Q_dc_S(:,2)./lstarSE; Q_bc_S(:,2)./lstarSE; Q_ac_S(:,2)./lstarSE; Q_am_S(:,2)./lstarSE];
% allz_S = [-10*RS/lstarSE; 10*RS/lstarSE; EARTH_dm_S(:,3)./lstarSE; MARS_dm_S(:,3)./lstarSE; Q_dmEM_S(:,3)./lstarSE; Q_dmSE_S(:,3)./lstarSE; Q_dc_S(:,3)./lstarSE; Q_bc_S(:,3)./lstarSE; Q_ac_S(:,3)./lstarSE; Q_am_S(:,3)./lstarSE];
% limx_S = axLimits(allx_S, pad);
% limy_S = axLimits(ally_S, pad);
% limz_S = axLimits(allz_S, pad);
% xlim(axS, limx_S)
% ylim(axS, limy_S)
% zlim(axS, limz_S)
% allx_EM = [(1-muEM-Rm/lstarEM).*lstarEM./lstarSE; (1-muEM+Rm/lstarEM).*lstarEM./lstarSE; transfer.departureOrbit.x.*lstarEM./lstarSE; q_dmEM_EM(:,1).*lstarEM./lstarSE];
% ally_EM = [-Rm/lstarSE; Rm/lstarSE; transfer.departureOrbit.y.*lstarEM./lstarSE; q_dmEM_EM(:,2).*lstarEM./lstarSE];
% allz_EM = [-Rm/lstarSE; Rm/lstarSE; transfer.departureOrbit.z.*lstarEM./lstarSE; q_dmEM_EM(:,3).*lstarEM./lstarSE];
% limx_EM = axLimits(allx_EM, pad);
% limy_EM = axLimits(ally_EM, pad);
% limz_EM = axLimits(allz_EM, pad);
% xlim(axEM, limx_EM)
% ylim(axEM, limy_EM)
% zlim(axEM, limz_EM)
% allx_SE = [1-muSE-RE/lstarSE; 1-muSE+RE/lstarSE; Moon_dm_SE(:,1); q_dmEM_SE(:,1); q_dmSE_SE(1:end-100,1)];
% ally_SE = [-RE/lstarSE; RE/lstarSE; Moon_dm_SE(:,2); q_dmEM_SE(:,2); q_dmSE_SE(1:end-100,2)];
% allz_SE = [-RE/lstarSE; RE/lstarSE; Moon_dm_SE(:,3); q_dmEM_SE(:,3); q_dmSE_SE(1:end-100,3)];
% limx_SE = axLimits(allx_SE, pad);
% limy_SE = axLimits(ally_SE, pad);
% limz_SE = axLimits(allz_SE, pad);
% xlim(axSE, limx_SE)
% ylim(axSE, limy_SE)
% zlim(axSE, limz_SE)
% shiftZLabel(axSE, hzSE, 0.15)
% allx_SM = [(1-muSM-10*RM/lstarSM).*lstarSM./lstarSE; (1-muSM+10*RM/lstarSM).*lstarSM./lstarSE; transfer.arrivalOrbit.x.*lstarSM./lstarSE; q_am_SM(101:end,1).*lstarSM./lstarSE];
% ally_SM = [-10*RM/lstarSE; 10*RM/lstarSE; transfer.arrivalOrbit.y.*lstarSM./lstarSE; q_am_SM(101:end,2).*lstarSM./lstarSE];
% allz_SM = [-10*RM/lstarSE; 10*RM/lstarSE; transfer.arrivalOrbit.z.*lstarSM./lstarSE; q_am_SM(101:end,3).*lstarSM./lstarSE];
% limx_SM = axLimits(allx_SM, pad);
% limy_SM = axLimits(ally_SM, pad);
% limz_SM = axLimits(allz_SM, pad);
% xlim(axSM, limx_SM)
% ylim(axSM, limy_SM)
% zlim(axSM, limz_SM)
% 
% % leg10 = legend(axEM, [Sun_S, Earth_EM, traj_Earth_S, Moon_EM, traj_Moon_SE, Venus_SV, traj_Venus_S, orbit_EM, orbit_SV, hist_dmEM_EM, hist_dc_S, hist_bc_S, hist_ac_S, hist_am_S, Deltav1], 'AutoUpdate', 'off');
% % drawnow;
% % set(leg10.EntryContainer.NodeChildren(end-1).Icon.Transform.Children.Children, 'ColorData', uint8([25; 25; 85; 255]));
% % legPos = [-0.25 0.2 0.3 0.7];
% % set(leg10, 'Units', 'normalized', 'Position', legPos);
% tit10 = annotation(fig10, 'textbox', [0.4 0.93 0.2 0.06], 'String', 'Elapsed Time: 0 d.', 'FontName', 'Times New Roman', 'FontSize', 18, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% 
% % v = VideoWriter('MarsMMATCR3BP_sample.mp4', 'MPEG-4');
% % v.FrameRate = 40;
% % v.Quality = 100;
% % open(v);
% 
% stride = 1;
% tailLen = Inf;
% setAxVisible(axS, true);
% setAxVisible(axEM, true);
% setAxVisible(axSE, false);
% setAxVisible(axSM, false);
% 
% N1 = length(t_dmEM_EM);
% for j = 1:stride:N1
%     i0 = max(1, j-tailLen);
% 
%     set(hist_dmEM_S, 'XData', Q_dmEM_S(i0:j,1)./lstarSE, 'YData', Q_dmEM_S(i0:j,2)./lstarSE, 'ZData', Q_dmEM_S(i0:j,3)./lstarSE);
%     set(mark_Earth_S, 'XData', sx_E+EARTH_dm_S(j,1)./lstarSE, 'YData', sy_E+EARTH_dm_S(j,2)./lstarSE, 'ZData', sz_E+EARTH_dm_S(j,3)./lstarSE);
%     set(mark_Mars_S, 'XData', sx_M+MARS_dm_S(j,1)./lstarSE, 'YData', sy_M+MARS_dm_S(j,2)./lstarSE, 'ZData', sz_M+MARS_dm_S(j,3)./lstarSE);
%     set(mark_dm_S, 'XData', Q_dmEM_S(j,1)./lstarSE, 'YData', Q_dmEM_S(j,2)./lstarSE, 'ZData', Q_dmEM_S(j,3)./lstarSE);
% 
%     set(hist_dmEM_EM, 'XData', q_dmEM_EM(i0:j,1).*lstarEM./lstarSE, 'YData', q_dmEM_EM(i0:j,2).*lstarEM./lstarSE, 'ZData', q_dmEM_EM(i0:j,3).*lstarEM./lstarSE);
%     set(mark_EM, 'XData', q_dmEM_EM(j,1).*lstarEM./lstarSE, 'YData', q_dmEM_EM(j,2).*lstarEM./lstarSE, 'ZData', q_dmEM_EM(j,3).*lstarEM./lstarSE);
% 
%     set(tit10, 'String', sprintf('Elapsed Time: %.0f d.', floor((t_all(j)-t_all(1))/3600/24)));
% 
%     if j == N1-1
%         setAxVisible(axEM, false);
%     end
%     
%     drawnow limitrate;
% %     set(leg10, 'Units', 'normalized', 'Position', legPos);
% 
%     pause(0.1)
% %     writeVideo(v, getframe(fig10));
% end
% stride = 2;
% set(hist_dmEM_SE, 'XData', q_dmEM_SE(:,1), 'YData', q_dmEM_SE(:,2), 'ZData', q_dmEM_SE(:,3));
% setAxVisible(axSE, true);
% N2 = length(t_dmSE_SE);
% for j = 2:stride:N2
%     i0 = max(1, j-tailLen);
% 
%     set(hist_dmSE_S, 'XData', Q_dmSE_S(i0:j,1)./lstarSE, 'YData', Q_dmSE_S(i0:j,2)./lstarSE, 'ZData', Q_dmSE_S(i0:j,3)./lstarSE);
%     set(mark_Earth_S, 'XData', sx_E+EARTH_dm_S(N1+j,1)./lstarSE, 'YData', sy_E+EARTH_dm_S(N1+j,2)./lstarSE, 'ZData', sz_E+EARTH_dm_S(N1+j,3)./lstarSE);
%     set(mark_Mars_S, 'XData', sx_M+MARS_dm_S(N1+j,1)./lstarSE, 'YData', sy_M+MARS_dm_S(N1+j,2)./lstarSE, 'ZData', sz_M+MARS_dm_S(N1+j,3)./lstarSE);
%     set(mark_dm_S, 'XData', Q_dmSE_S(j,1)./lstarSE, 'YData', Q_dmSE_S(j,2)./lstarSE, 'ZData', Q_dmSE_S(j,3)./lstarSE);
% 
%     set(hist_dmSE_SE, 'XData', q_dmSE_SE(i0:j,1), 'YData', q_dmSE_SE(i0:j,2), 'ZData', q_dmSE_SE(i0:j,3));
%     set(mark_Moon_SE, 'XData', sx_m+Moon_dm_SE(N1+j,1), 'YData', sy_m+Moon_dm_SE(N1+j,2), 'ZData', sz_m+Moon_dm_SE(N1+j,3));
%     set(mark_SE, 'XData', q_dmSE_SE(j,1), 'YData', q_dmSE_SE(j,2), 'ZData', q_dmSE_SE(j,3));
% 
%     set(tit10, 'String', sprintf('Elapsed Time: %.0f d.', floor((t_all(N1+j)-t_all(1))/3600/24)));
% 
%     if j >= N2-95
%         setAxVisible(axSE, false);
%     end
% 
%     drawnow limitrate;
% %     set(leg10, 'Units', 'normalized', 'Position', legPos);
% 
%     pause(0.1)
% %     writeVideo(v, getframe(fig10));
% end
% N3 = length(t_dc_S);
% for j = 2:stride:N3
%     i0 = max(1, j-tailLen);
% 
%     set(hist_dc_S, 'XData', Q_dc_S(i0:j,1)./lstarSE, 'YData', Q_dc_S(i0:j,2)./lstarSE, 'ZData', Q_dc_S(i0:j,3)./lstarSE);
%     set(mark_Earth_S, 'XData', sx_E+EARTH_dm_S(N1+N2+j,1)./lstarSE, 'YData', sy_E+EARTH_dm_S(N1+N2+j,2)./lstarSE, 'ZData', sz_E+EARTH_dm_S(N1+N2+j,3)./lstarSE);
%     set(mark_Mars_S, 'XData', sx_M+MARS_dm_S(N1+N2+j,1)./lstarSE, 'YData', sy_M+MARS_dm_S(N1+N2+j,2)./lstarSE, 'ZData', sz_M+MARS_dm_S(N1+N2+j,3)./lstarSE);
%     set(mark_dm_S, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
%     set(mark_dc_S, 'XData', Q_dc_S(j,1)./lstarSE, 'YData', Q_dc_S(j,2)./lstarSE, 'ZData', Q_dc_S(j,3)./lstarSE);
% 
%     set(tit10, 'String', sprintf('Elapsed Time: %.0f d.', floor((t_all(N1+N2+j)-t_all(1))/3600/24)));
% 
%     drawnow limitrate;
% %     set(leg10, 'Units', 'normalized', 'Position', legPos);
% 
%     pause(0.1)
% %     writeVideo(v, getframe(fig10));
% end
% N4 = length(t_bc_S);
% set(Deltav1, 'XData', Q_bc_S(1,1)./lstarSE, 'YData', Q_bc_S(1,2)./lstarSE, 'ZData', Q_bc_S(1,3)./lstarSE);
% for j = 2:stride:N4
%     i0 = max(1, j-tailLen);
% 
%     set(hist_bc_S, 'XData', Q_bc_S(i0:j,1)./lstarSE, 'YData', Q_bc_S(i0:j,2)./lstarSE, 'ZData', Q_bc_S(i0:j,3)./lstarSE);
%     set(mark_Earth_S, 'XData', sx_E+EARTH_dm_S(N1+N2+N3+j,1)./lstarSE, 'YData', sy_E+EARTH_dm_S(N1+N2+N3+j,2)./lstarSE, 'ZData', sz_E+EARTH_dm_S(N1+N2+N3+j,3)./lstarSE);
%     set(mark_Mars_S, 'XData', sx_M+MARS_dm_S(N1+N2+N3+j,1)./lstarSE, 'YData', sy_M+MARS_dm_S(N1+N2+N3+j,2)./lstarSE, 'ZData', sz_M+MARS_dm_S(N1+N2+N3+j,3)./lstarSE);
%     set(mark_dc_S, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
%     set(mark_bc_S, 'XData', Q_bc_S(j,1)./lstarSE, 'YData', Q_bc_S(j,2)./lstarSE, 'ZData', Q_bc_S(j,3)./lstarSE);
% 
%     set(tit10, 'String', sprintf('Elapsed Time: %.0f d.', floor((t_all(N1+N2+N3+j)-t_all(1))/3600/24)));
% 
%     drawnow limitrate;
% %     set(leg10, 'Units', 'normalized', 'Position', legPos);
% 
%     pause(0.1)
% %     writeVideo(v, getframe(fig10));
% end
% N5 = length(t_ac_S);
% set(Deltav2, 'XData', Q_ac_S(1,1)./lstarSE, 'YData', Q_ac_S(1,2)./lstarSE, 'ZData', Q_ac_S(1,3)./lstarSE);
% for j = 2:stride:N5
%     i0 = max(1, j-tailLen);
% 
%     set(hist_ac_S, 'XData', Q_ac_S(i0:j,1)./lstarSE, 'YData', Q_ac_S(i0:j,2)./lstarSE, 'ZData', Q_ac_S(i0:j,3)./lstarSE);
%     set(mark_Earth_S, 'XData', sx_E+EARTH_dm_S(N1+N2+N3+N4+j,1)./lstarSE, 'YData', sy_E+EARTH_dm_S(N1+N2+N3+N4+j,2)./lstarSE, 'ZData', sz_E+EARTH_dm_S(N1+N2+N3+N4+j,3)./lstarSE);
%     set(mark_Mars_S, 'XData', sx_M+MARS_dm_S(N1+N2+N3+N4+j,1)./lstarSE, 'YData', sy_M+MARS_dm_S(N1+N2+N3+N4+j,2)./lstarSE, 'ZData', sz_M+MARS_dm_S(N1+N2+N3+N4+j,3)./lstarSE);
%     set(mark_bc_S, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
%     set(mark_ac_S, 'XData', Q_ac_S(j,1)./lstarSE, 'YData', Q_ac_S(j,2)./lstarSE, 'ZData', Q_ac_S(j,3)./lstarSE);
% 
%     set(tit10, 'String', sprintf('Elapsed Time: %.0f d.', floor((t_all(N1+N2+N3+N4+j)-t_all(1))/3600/24)));
% 
%     drawnow limitrate;
% %     set(leg10, 'Units', 'normalized', 'Position', legPos);
% 
%     pause(0.1)
% %     writeVideo(v, getframe(fig10));
% end
% N6 = length(t_am_SM);
% for j = 2:stride:N6
%     i0 = max(1, j-tailLen);
% 
%     set(hist_am_S, 'XData', Q_am_S(i0:j,1)./lstarSE, 'YData', Q_am_S(i0:j,2)./lstarSE, 'ZData', Q_am_S(i0:j,3)./lstarSE);
%     set(mark_Earth_S, 'XData', sx_E+EARTH_dm_S(N1+N2+N3+N4+N5+j,1)./lstarSE, 'YData', sy_E+EARTH_dm_S(N1+N2+N3+N4+N5+j,2)./lstarSE, 'ZData', sz_E+EARTH_dm_S(N1+N2+N3+N4+N5+j,3)./lstarSE);
%     set(mark_Mars_S, 'XData', sx_M+MARS_dm_S(N1+N2+N3+N4+N5+j,1)./lstarSE, 'YData', sy_M+MARS_dm_S(N1+N2+N3+N4+N5+j,2)./lstarSE, 'ZData', sz_M+MARS_dm_S(N1+N2+N3+N4+N5+j,3)./lstarSE);
%     set(mark_ac_S, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
%     set(mark_am_S, 'XData', Q_am_S(j,1)./lstarSE, 'YData', Q_am_S(j,2)./lstarSE, 'ZData', Q_am_S(j,3)./lstarSE);
% 
%     set(hist_amSM_SM, 'XData', q_am_SM(i0:j,1).*lstarSM./lstarSE, 'YData', q_am_SM(i0:j,2).*lstarSM./lstarSE, 'ZData', q_am_SM(i0:j,3).*lstarSM./lstarSE);
%     set(mark_SM, 'XData', q_am_SM(j,1).*lstarSM./lstarSE, 'YData', q_am_SM(j,2).*lstarSM./lstarSE, 'ZData', q_am_SM(j,3).*lstarSM./lstarSE);
% 
%     set(tit10, 'String', sprintf('Elapsed Time: %.0f d.', floor((t_all(N1+N2+N3+N4+N5+j)-t_all(1))/3600/24)));
% 
%     if j >= 95
%         setAxVisible(axSM, true);
%     end
% 
%     drawnow limitrate;
% %     set(leg10, 'Units', 'normalized', 'Position', legPos);
% 
%     pause(0.1)
% %     writeVideo(v, getframe(fig10));
% end
% 
% % close(v);

%% Functions
function setAxVisible(ax, tf)
    vis = matlab.lang.OnOffSwitchState(tf);
    ax.Visible = vis;
    set(findall(ax), 'Visible', vis);
end