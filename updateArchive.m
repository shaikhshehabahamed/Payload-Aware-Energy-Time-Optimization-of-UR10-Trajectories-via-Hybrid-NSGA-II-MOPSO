function [archivePos, archiveCosts] = updateArchive( ...
    archivePos, archiveCosts, newPos, newCosts, archiveMax)
%UPDATEARCHIVE Update external archive with nondominated solutions.

% Concatenate
allPos   = [archivePos; newPos];
allCosts = [archiveCosts; newCosts];

% Keep only finite-cost solutions
validIdx = all(isfinite(allCosts),2);
allPos   = allPos(validIdx,:);
allCosts = allCosts(validIdx,:);

if isempty(allPos)
    archivePos   = [];
    archiveCosts = [];
    return;
end

% Non-dominated sorting on all solutions
fronts = nondominatedSort(allCosts);
F1 = fronts{1};
ndPos   = allPos(F1,:);
ndCosts = allCosts(F1,:);

% If archive too large, use crowding distance to truncate
if size(ndPos,1) > archiveMax
    cd = computeCrowdingDistance(ndCosts);
    [~, sortIdx] = sort(cd, 'descend');
    keepIdx = sortIdx(1:archiveMax);
    ndPos   = ndPos(keepIdx,:);
    ndCosts = ndCosts(keepIdx,:);
end

archivePos   = ndPos;
archiveCosts = ndCosts;

end