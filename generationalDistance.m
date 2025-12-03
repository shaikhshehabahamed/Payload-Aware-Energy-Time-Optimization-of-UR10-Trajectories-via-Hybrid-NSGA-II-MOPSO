function gd = generationalDistance(front, refFront)
%GENERATIONALDISTANCE Generational distance between front and reference.
%   front: N x M
%   refFront: K x M (approximate true Pareto front)

if isempty(front) || isempty(refFront)
    gd = Inf;
    return;
end

N = size(front,1);
sumSq = 0;

for i = 1:N
    diffs = refFront - front(i,:);
    dists = sqrt(sum(diffs.^2, 2));
    dmin = min(dists);
    sumSq = sumSq + dmin^2;
end

gd = sqrt(sumSq)/N;

end