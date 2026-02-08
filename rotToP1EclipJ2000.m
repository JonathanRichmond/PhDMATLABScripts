%%% rotToP1EclipJ2000
%%% Jonathan LeFevre Richmond
%%% C: 29 May 2024
%%% U: 7 February 2026

function P1EclipJ2000States = rotToP1EclipJ2000(mu, initialEpoch, P1, gm1, P2, lstar, tstar, times, states)

cspice_furnsh({'naif0012.tls', 'de430.bsp', 'de440.bsp', 'mar097.bsp'}); % Load kernels

[P2InitialStateDim, ~] = getSPICEEclipJ2000(initialEpoch, 0, P2, P1); % Body initial state vector [dim]
initialEpochTime = cspice_str2et(initialEpoch);
P2SPICEElements = cspice_oscltx(P2InitialStateDim, initialEpochTime, gm1)'; % Body initial orbital elements [dim]
if P2(1) == 'E'
    P2SPICEElements(3) = 0;
end
timesDim = times.*tstar;

P1EclipJ2000States = zeros(length(times), 6);
for i = 1:length(times)
    stateDim = [states(i,1:3).*lstar, states(i,4:6).*lstar./tstar];
    statePrimaryDim = stateDim-[-mu*lstar, 0, 0, 0, 0, 0];
    P2Elements = [lstar, 0, P2SPICEElements(3:5), P2SPICEElements(9)+times(i), initialEpochTime+timesDim(i), P2SPICEElements(8)]';
    P2StateDim = cspice_conics(P2Elements, initialEpochTime+timesDim(i));
    xhat = P2StateDim(1:3)./lstar;
    zhat = cross(P2StateDim(1:3), P2StateDim(4:6))./norm(cross(P2StateDim(1:3), P2StateDim(4:6)));
    yhat = cross(zhat, xhat);
    C = [xhat, yhat, zhat];
    thetadotDim = 1/tstar;
    Cdot = [thetadotDim.*yhat, -1*thetadotDim.*xhat, zeros(3, 1)];
    N = [C, zeros(3, 3); Cdot C];
    P1EclipJ2000StateDim = (N*statePrimaryDim')';
    P1EclipJ2000States(i,:) = [P1EclipJ2000StateDim(1:3)./lstar, P1EclipJ2000StateDim(4:6).*tstar./lstar];
end

cspice_kclear; % Clear kernels