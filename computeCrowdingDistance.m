function cd = computeCrowdingDistance(costs)
%COMPUTECROWDINGDISTANCE Crowding distance for a single front

[n, nObj] = size(costs);
cd = zeros(n,1);

if n == 0
    return;
elseif n == 1
    cd(1) = Inf;
    return;
elseif n == 2
    cd(:) = Inf;
    return;
end

for j = 1:nObj
    [sortedVals, sortIdx] = sort(costs(:,j));
    cd(sortIdx(1))   = Inf;
    cd(sortIdx(end)) = Inf;

    range = sortedVals(end) - sortedVals(1);
    if range == 0
        range = 1; % avoid division by zero
    end

    for i = 2:(n-1)
        cd(sortIdx(i)) = cd(sortIdx(i)) + ...
            (sortedVals(i+1) - sortedVals(i-1))/range;
    end
end

end