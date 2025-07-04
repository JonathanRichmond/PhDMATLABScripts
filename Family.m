%%% Plot Family
%%% Jonathan Richmond
%%% C: 3 July 2025

clear
load('../PhDScripts/FamilyData/CR3BPEMLFRPros.mat')
k = length(LFRPros.trajs);

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

%% Libration Points
g1 = mu;
delg1 = 1;
while abs(delg1) > eps
    f = ((1-mu)/((1-g1)^2))-(mu/(g1^2))-1+mu+g1;
    fprime = ((2*(1-mu))/((1-g1)^3))+((2*mu)/(g1^3))+1;
    g1New = g1-(f/fprime);
    delg1 = g1New-g1;
    g1 = g1New;
end
a1 = 1-mu-g1;

g2 = mu;
delg2 = 1;
while abs(delg2) > eps
    f = ((1-mu)/((1+g2)^2))+(mu/(g2^2))-1+mu-g2;
    fprime = ((-2*(1-mu))/((1+g2)^3))-((2*mu)/(g2^3))-1;
    g2New = g2-(f/fprime);
    delg2 = g2New-g2;
    g2 = g2New;
end
a2 = 1-mu+g2;

g3 = mu;
delg3 = 1;
while abs(delg3) > eps
    f = ((1-mu)/(g3^2))+(mu/((1+g3)^2))-mu-g3;
    fprime = ((-2*(1-mu))/(g3^3))-((2*mu)/((1+g3)^3))-1;
    g3New = g3-(f/fprime);
    delg3 = g3New-g3;
    g3 = g3New;
end
a3 = -1*mu-g3;

a45 = 0.5-mu;
b4 = sqrt(3)/2;
b5 = -b4;

%% Plot
fig1 = figure("Position", [200 100 1200 750]);
hold on
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
scatter3(a1, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "$L_{1}$")
scatter3(a2, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "$L_{2}$")
% scatter3(a3, 0, 0, 20, 'g', 'filled', 'd', 'DisplayName', "$L_{3}$")
scatter3(a45, b4, 0, 20, 'b', 'filled', 'd', 'DisplayName', "$L_{4}$")
% scatter3(a45, b5, 0, 20, [1 0 1], 'filled', 'd', 'DisplayName', "$L_{5}$")
for j = 1:40:k
    traj = plot3(LFRPros.trajs{j}.x, LFRPros.trajs{j}.y, LFRPros.trajs{j}.z, 'b', 'HandleVisibility', 'off', 'DisplayName', "Lunar Free Return");
    if j == 1
        set(traj, 'HandleVisibility', 'on')
    end
end
axis equal
% axis
axis([-0.12 1.24 -0.05 0.87])
grid on
xlabel("$x$ [EM ndim]", 'Interpreter', 'latex')
ylabel("$y$ [EM ndim]", 'Interpreter', 'latex')
zlabel("$z$ [EM ndim]", 'Interpreter', 'latex')
title("Earth-Moon Rot.", 'Interpreter', 'latex')
legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'Color', 'k');
view(2)
hold off
% exportgraphics(fig1, 'Family_1.png', 'BackgroundColor', 'k')