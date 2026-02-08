%%% planetRotToSunEclipJ2000
%%% Jonathan LeFevre Richmond
%%% C: 30 October 2023
%%% U: 7 February 2026

function SunEclipJ2000States = planetRotToSunEclipJ2000(mu, initialEpoch, body, lstar, tstar, times, states)

cspice_furnsh({'naif0012.tls', 'de430.bsp', 'de440.bsp', 'mar097.bsp'}); % Load kernels
gmS = 1.3271244004193930E11; % Sun gravitational parameter [km^3/s^2]

[bodyInitialStateDim, ~] = getSPICEEclipJ2000(initialEpoch, 0, body, 'Sun'); % Body initial state vector [dim]
initialEpochTime = cspice_str2et(initialEpoch);
bodySPICEElements = cspice_oscltx(bodyInitialStateDim, initialEpochTime, gmS)'; % Body initial orbital elements [dim]
if body(1) == 'E'
    bodySPICEElements(3) = 0;
end
timesDim = times.*tstar;

SunEclipJ2000States = zeros(length(times), 6);
for i = 1:length(times)
    stateDim = [states(i,1:3).*lstar, states(i,4:6).*lstar./tstar];
    statePrimaryDim = stateDim-[(-1*mu)*lstar, 0, 0, 0, 0, 0];
    bodyElements = [lstar, 0, bodySPICEElements(3:5), bodySPICEElements(9)+times(i), initialEpochTime+timesDim(i), bodySPICEElements(8)]';
    bodyStateDim = cspice_conics(bodyElements, initialEpochTime+timesDim(i));
    xhat = bodyStateDim(1:3)./lstar;
    zhat = cross(bodyStateDim(1:3), bodyStateDim(4:6))./norm(cross(bodyStateDim(1:3), bodyStateDim(4:6)));
    yhat = cross(zhat, xhat);
    C = [xhat, yhat, zhat];
    thetadotDim = 1/tstar;
    Cdot = [thetadotDim.*yhat, -1*thetadotDim.*xhat, zeros(3, 1)];
    N = [C, zeros(3, 3); Cdot C];
    SunEclipJ2000StateDim = (N*statePrimaryDim')';
    SunEclipJ2000States(i,:) = [SunEclipJ2000StateDim(1:3)./lstar, SunEclipJ2000StateDim(4:6).*tstar./lstar];
end

cspice_kclear; % Clear kernels
