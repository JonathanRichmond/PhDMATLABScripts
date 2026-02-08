%%% P2EclipJ2000ToRot
%%% Jonathan LeFevre Richmond
%%% C: 29 May 2024
%%% U: 7 February 2026

function RotStates = P2EclipJ2000ToRot(mu, initialEpoch, P1, gm1, P2, lstar, tstar, times, states)

cspice_furnsh({'naif0012.tls', 'de430.bsp', 'de440.bsp', 'mar097.bsp'}); % Load kernels

[P2InitialStateDim, ~] = getSPICEEclipJ2000(initialEpoch, 0, P2, P1); % Body initial state vector [dim]
initialEpochTime = cspice_str2et(initialEpoch);
P2SPICEElements = cspice_oscltx(P2InitialStateDim, initialEpochTime, gm1)'; % Body initial orbital elements [dim]
if P2(1) == 'E'
    P2SPICEElements(3) = 0;
end
timesDim = times.*tstar;

RotStates = zeros(length(times), 6);
for i = 1:length(times)
    P2EclipJ2000StateDim = [states(i,1:3).*lstar, states(i,4:6).*lstar./tstar];
    P2Elements = [lstar, 0, P2SPICEElements(3:5), P2SPICEElements(9)+times(i), initialEpochTime+timesDim(i), P2SPICEElements(8)]';
    P2StateDim = cspice_conics(P2Elements, initialEpochTime+timesDim(i));
    P1EclipJ2000StateDim = P2StateDim'+P2EclipJ2000StateDim;
    xhat = P2StateDim(1:3)./lstar;
    zhat = cross(P2StateDim(1:3), P2StateDim(4:6))./norm(cross(P2StateDim(1:3), P2StateDim(4:6)));
    yhat = cross(zhat, xhat);
    C = [xhat, yhat, zhat];
    thetadotDim = 1/tstar;
    Cdot = [thetadotDim.*yhat, -1*thetadotDim.*xhat, zeros(3, 1)];
    N = [C, zeros(3, 3); Cdot C];
    P1RotStateDim = (N\P1EclipJ2000StateDim')';
    P1RotState = [P1RotStateDim(1:3)./lstar, P1RotStateDim(4:6).*tstar./lstar];
    RotStates(i,:) = P1RotState+[-mu, 0, 0, 0, 0, 0];
end

cspice_kclear; % Clear kernels