function hv = hypervolumeMetric(front, refPoint)
%HYPERVOLUMEMETRIC Hypervolume in 2D (minimization).
%   front: N x 2 objective matrix
%   refPoint: 1 x 2 reference point (worse than all fronts)

if isempty(front)
    hv = 0;
    return;
end

% Remove non-finite
front = front(all(isfinite(front),2),:);
if isempty(front)
    hv = 0;
    return;
end

% Sort by first objective (e.g., energy)
frontSorted = sortrows(front,1);

% Ensure monotonic non-dominated set (f2 decreasing)
pruned = [];
bestF2 = Inf;
for i = 1:size(frontSorted,1)
    f = frontSorted(i,:);
    if f(2) < bestF2
        pruned = [pruned; f]; %#ok<AGROW>
        bestF2 = f(2);
    end
end

front = pruned;
hv = 0;
prevF2 = refPoint(2);

for i = 1:size(front,1)
    width  = max(refPoint(1) - front(i,1), 0);
    height = max(prevF2 - front(i,2), 0);
    hv = hv + width*height;
    prevF2 = front(i,2);
end

end