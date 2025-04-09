%%% getHEM
%%% Jonathan Richmond
%%% C: 8 April 2025

function HEM = getHEM(muEM, mS, aS, q)

v2 = q(4)^2+q(5)^2+q(6)^2;
r_13 = sqrt((q(1)+muEM)^2+q(2)^2+q(3)^2);
r_23 = sqrt((q(1)-1+muEM)^2+q(2)^2+q(3)^2);
r_43 = sqrt((q(1)-aS*cos(q(7)))^2+(q(2)-aS*sin(q(7)))^2+q(3)^2);
U = (1-muEM)/r_13+muEM/r_23+mS/r_43-mS*(q(1)*cos(q(7))+q(2)*sin(q(7)))/aS^2+(1/2)*(q(1)^2+q(2)^2);

HEM = 2*U-v2;