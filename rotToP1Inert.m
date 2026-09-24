%%% rotToP1Inert
%%% Jonathan LeFevre Richmond
%%% C: 23 September 2026

function P1InertStates = rotToP1Inert(mu, times, states)

statesPrimary = states-[-mu, 0, 0, 0, 0, 0];
P1InertStates = zeros(length(times), 6);
for j = 1:length(times)
    t = times(j);
    C = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
    Cdot = [-sin(t), -cos(t), 0; cos(t), -sin(t), 0; 0, 0, 0];
    N = [C, zeros(3, 3); Cdot C];
    P1InertStates(j,:) = (N*statesPrimary(j,:)')';
end