function results = hybridOptimization(problem, mopsoParams, nsga2Params)
%HYBRIDOPTIMIZATION Run MOPSO, NSGA-II and Hybrid NSGA-II+MOPSO
%   results.MOPSO, results.NSGA2, results.Hybrid
%   results.metrics.(Algorithm).HV, GD, Spread

fprintf('  -> Running MOPSO only...\n');
mopsoRes = mopsoPhase(problem, mopsoParams);

fprintf('  -> Running NSGA-II only...\n');
nsga2Res = nsga2Phase(problem, nsga2Params, []);

fprintf('  -> Running Hybrid (MOPSO -> NSGA-II)...\n');
% Seed NSGA-II with MOPSO archive solutions
initialPositionsHybrid = mopsoRes.archivePositions;
hybridRes = nsga2Phase(problem, nsga2Params, initialPositionsHybrid);

% Build reference front from union of all algorithms
allFronts = [mopsoRes.paretoFront; nsga2Res.paretoFront; hybridRes.paretoFront];
validIdx  = all(all(isfinite(allFronts)),2);
allFronts = allFronts(validIdx,:);

if isempty(allFronts)
    warning('No feasible solutions found for any algorithm.');
    refFront = [];
    refPoint = [1,1]; % dummy
else
    % Reference front = nondominated front of union
    refFrontIdx = nondominatedSort(allFronts);
    refFront    = allFronts(refFrontIdx{1},:);

    % Reference point for HV: slightly worse than worst in union
    maxVals = max(allFronts,[],1);
    refPoint = 1.1*maxVals;
end

% Compute metrics for each algorithm
metrics.MOPSO = computeMetricsForAlgo(mopsoRes.paretoFront, refFront, refPoint);
metrics.NSGA2 = computeMetricsForAlgo(nsga2Res.paretoFront, refFront, refPoint);
metrics.Hybrid= computeMetricsForAlgo(hybridRes.paretoFront, refFront, refPoint);

% Package results
results.MOPSO  = mopsoRes;
results.NSGA2  = nsga2Res;
results.Hybrid = hybridRes;
results.metrics= metrics;

end

function m = computeMetricsForAlgo(front, refFront, refPoint)

if isempty(front) || isempty(refFront)
    m.HV     = 0;
    m.GD     = Inf;
    m.Spread = Inf;
    return;
end

m.HV     = hypervolumeMetric(front, refPoint);
m.GD     = generationalDistance(front, refFront);
m.Spread = spreadMetric(front, refFront);

end