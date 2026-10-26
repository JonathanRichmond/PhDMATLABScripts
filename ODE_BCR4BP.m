%%% ODE_BCR4BP
%%% Jonathan LeFevre Richmond
%%% C: 14 October 2026
%%% U: 15 October 2026

function dqdt = ODE_BCR4BP(~, q, mu, m4, a4, mS, mstar)

x = q(1);
y = q(2);
z = q(3);
xdot = q(4);
ydot = q(5);
zdot = q(6);
thetaS = q(7);

r13 = sqrt((x+mu)^2+y^2+z^2);
r23 = sqrt((x-1+mu)^2+y^2+z^2);
r43 = sqrt((x-a4*cos(thetaS))^2+(y-a4*sin(thetaS))^2+z^2);

xddot = 2*ydot+x-(1-mu)*(x+mu)/(r13^3)-mu*(x-1+mu)/(r23^3)-m4*(x-a4*cos(thetaS))/(r43^3)-m4*cos(thetaS)/(a4^2);
yddot = -2*xdot+y-(1-mu)*y/(r13^3)-mu*y/(r23^3)-m4*(y-a4*sin(thetaS))/(r43^3)-m4*sin(thetaS)/(a4^2);
zddot = -1*(1-mu)*z/(r13^3)-mu*z/(r23^3)-m4*z/(r43^3);

thetaSdot = sqrt((mS/mstar+1)/(a4^3))-1;

dqdt = [xdot; ydot; zdot; xddot; yddot; zddot; thetaSdot];