%%% ODE_C
%%% Jonathan LeFevre Richmond
%%% C: 20 December 2023
%%% U: 31 January 2026

function dqdt = ODE_C(~, y, mu)

r = y(1:3);
v = y(4:6);
rMag = norm(r);

dqdt = [v; -(mu/rMag^3)*r];