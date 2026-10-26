%%% ODE_CR3BP_JC
%%% Jonathan LeFevre Richmond
%%% C: 22 October 2026

function dqdt = ODE_CR3BP_JC(~, q, mu, m4, a4, mS, mstar)

x = q(1);
y = q(2);
z = q(3);
xdot = q(4);
ydot = q(5);
zdot = q(6);
thetaS = q(7);

r13 = sqrt((x+mu)^2+y^2+z^2);
r23 = sqrt((x-1+mu)^2+y^2+z^2);

xddot = 2*ydot+x-(1-mu)*(x+mu)/(r13^3)-mu*(x-1+mu)/(r23^3);
yddot = -2*xdot+y-(1-mu)*y/(r13^3)-mu*y/(r23^3);
zddot = -1*(1-mu)*z/(r13^3)-mu*z/(r23^3);

thetaSdot = sqrt((mS/mstar+1)/(a4^3))-1;
JCdot = -(m4/(a4^3))*(x*xdot*(1+3*cos(2*thetaS))+3*y*xdot*sin(2*thetaS)+3*x*ydot*sin(2*thetaS)+y*ydot*(1-3*cos(2*thetaS)));

dqdt = [xdot; ydot; zdot; xddot; yddot; zddot; thetaSdot; JCdot];