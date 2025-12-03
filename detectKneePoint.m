function [kneeIndex, kneePoint] = detectKneePoint(paretoFront)
%DETECTKNEEPOINT Simple knee-point detection via max distance to chord.
%   paretoFront: N x 2 matrix [energy, time] (both minimized)

if isempty(paretoFront)
    kneeIndex = [];
    kneePoint = [];
    return;
end

costs = paretoFront;
minVals = min(costs,[],1);
maxVals = max(costs,[],1);
normCosts = (costs - minVals)./(maxVals - minVals + eps);

% Choose two extreme points: min f1 and min f2
[~, idx1] = min(normCosts(:,1));
[~, idx2] = min(normCosts(:,2));
p1 = normCosts(idx1,:);
p2 = normCosts(idx2,:);

if all(abs(p1 - p2) < 1e-12)
    % Degenerate case: choose point with minimal sum of norms
    [~, kneeIndex] = min(sum(normCosts,2));
    kneePoint = paretoFront(kneeIndex,:);
    return;
end

v = p2 - p1;
vNorm = sqrt(sum(v.^2));

N = size(normCosts,1);
distances = zeros(N,1);

for i = 1:N
    w = normCosts(i,:) - p1;
    % Perpendicular distance in 2D
    distances(i) = abs(v(1)*w(2) - v(2)*w(1))/vNorm;
end

[~, kneeIndex] = max(distances);
kneePoint = paretoFront(kneeIndex,:);

end