%%% getSPICEEclipJ2000
%%% Jonathan LeFevre Richmond
%%% C: 30 October 2023
%%% U: 18 March 2025

function [state, ltime] = getSPICEEclipJ2000(startstr, tdim, targetstr, obsstr)

start = cspice_str2et(startstr); % Start epoch (Julian)
et = tdim'+start; % Epochs (Julian)

[state, ltime] = cspice_spkezr(targetstr, et, 'ECLIPJ2000', 'NONE', obsstr); % State vectors