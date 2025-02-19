%%% BCR4BPDev.jl
%%% Jonathan Richmond
%%% C: 19 February 2025

clear
load("../PhDScripts/Output/BCR4BPDev.mat")

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
g1(1) = mu;
delg1 = 1;
i = 0;
while abs(delg1) > eps
    f = ((1-mu)/((1-g1(i+1))^2))-(mu/(g1(i+1)^2))-1+mu+g1(i+1);
    fprime = ((2*(1-mu))/((1-g1(i+1))^3))+((2*mu)/(g1(i+1)^3))+1;
    g1(i+2) = g1(i+1)-(f/fprime);
    delg1 = g1(i+2)-g1(i+1);
    i = i+1;
end
a1 = 1-mu-g1(end);

g2(1) = mu;
delg2 = 1;
i = 0;
while abs(delg2) > eps
    f = ((1-mu)/((1+g2(i+1))^2))+(mu/(g2(i+1)^2))-1+mu-g2(i+1);
    fprime = ((-2*(1-mu))/((1+g2(i+1))^3))-((2*mu)/(g2(i+1)^3))-1;
    g2(i+2) = g2(i+1)-(f/fprime);
    delg2 = g2(i+2)-g2(i+1);
    i = i+1;
end
a2 = 1-mu+g2(end);

g3(1) = mu;
delg3 = 1;
i = 0;
while abs(delg3) > eps
    f = ((1-mu)/(g3(i+1)^2))+(mu/((1+g3(i+1))^2))-mu-g3(i+1);
    fprime = ((-2*(1-mu))/(g3(i+1)^3))-((2*mu)/((1+g3(i+1))^3))-1;
    g3(i+2) = g3(i+1)-(f/fprime);
    delg3 = g3(i+2)-g3(i+1);
    i = i+1;
end
a3 = -1*mu-g3(end);

a45 = 0.5-mu;
b4 = sqrt(3)/2;
b5 = -b4;

%% Plot
fig1 = figure("Position", [200 100 1200 750]);
Earth = plot3DBody("Earth", RE/lstar, [-mu, 0, 0]);
set(Earth, 'DisplayName', "Earth")
hold on
Moon = plot3DBody("Moon", Rm/lstar, [1-mu, 0, 0]);
set(Moon, 'DisplayName', "Moon")
scatter3(a1, 0, 0, 20, 'r', 'filled', 'd', 'DisplayName', "L_{1}")
scatter3(a2, 0, 0, 20, [1 0.5 0], 'filled', 'd', 'DisplayName', "L_{2}")
scatter3(a3, 0, 0, 20, 'g', 'filled', 'd', 'DisplayName', "L_{3}")
scatter3(a45, b4, 0, 20, 'b', 'filled', 'd', 'DisplayName', "L_{4}")
scatter3(a45, b5, 0, 20, [1 0 1], 'filled', 'd', 'DisplayName', "L_{5}")
plot3(trajCR3BP.x, trajCR3BP.y, trajCR3BP.z, 'DisplayName', "CR3BP Orbit")
scatter3(trajBCR4BP12.x, trajBCR4BP12.y, trajBCR4BP12.z, 10*ones(length(trajBCR4BP12.t), 1), SunAngleColor(trajBCR4BP12.theta4), 'filled', 'DisplayName', 'BCR4BP E-M Traj.')
axis equal
grid on
xlabel("x [EM ndim]")
ylabel("y [EM ndim]")
zlabel("z [EM ndim]")
title("L_{1} Lyapunov: \theta_{S} = "+num2str(rad2deg(trajBCR4BP12.theta4(1)), '%.0f')+" deg")
leg = legend('Location', 'northeastoutside');
phasemap;
phasebar('Location', 'northeast', 'Size', 0.25)
set(gca, 'Color', 'k');
hold off
exportgraphics(fig1, 'BCR4BPDev_1.png', 'BackgroundColor', 'k')