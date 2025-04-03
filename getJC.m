%%% getJC
%%% Jonathan Richmond
%%% C: 2 April 2025

function JC = getJC(mu, q)

v2 = q(4)^2+q(5)^2+q(6)^2;
r_13 = sqrt((q(1)+mu)^2+q(2)^2+q(3)^2);
r_23 = sqrt((q(1)-1+mu)^2+q(2)^2+q(3)^2);
U = (1-mu)/r_13+mu/r_23+(1/2)*(q(1)^2+q(2)^2);

JC = 2*U-v2;