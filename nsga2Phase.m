function nsgaRes = nsga2Phase(problem, params, initialPositions)
%NSGA2PHASE NSGA-II refinement phase.
%   nsgaRes.population  : final population struct
%   nsgaRes.paretoFront : N x M matrix of objective values
%   nsgaRes.paretoSet   : N x nVar matrix of decision vectors
%
%   initialPositions: optional matrix [nInit x nVar] used for seeding
%   the initial population with solutions from MOPSO etc.

nVar   = problem.nVar;
varMin = problem.varMin;
varMax = problem.varMax;

nPop   = params.nPop;
maxGen = params.maxGen;
pC     = params.pCrossover;
pM     = params.pMutation;
etaC   = params.etaC;
etaM   = params.etaM;

% Individual template
empty_ind.Position = [];
empty_ind.Cost     = [];
empty_ind.Rank     = [];
empty_ind.CrowdingDistance = [];

pop = repmat(empty_ind, nPop, 1);

% ----- Initialize population -----
if nargin < 3 || isempty(initialPositions)
    nInit = 0;
else
    nInit = min(size(initialPositions,1), nPop);
end

% Seed with initial positions (if provided)
for i = 1:nInit
    pop(i).Position = initialPositions(i,:);
    pop(i).Cost     = objectiveFunction(pop(i).Position, problem);
end

% Fill remaining with random individuals
for i = (nInit+1):nPop
    pop(i).Position = varMin + rand(1,nVar).*(varMax - varMin);
    pop(i).Cost     = objectiveFunction(pop(i).Position, problem);
end

% ----- Main NSGA-II loop -----
for gen = 1:maxGen

    % --- Non-dominated sorting ---
    costMatrix = vertcat(pop.Cost);
    fronts = nondominatedSort(costMatrix);

    % --- Crowding distance & rank assignment ---
    for f = 1:numel(fronts)
        frontIdx = fronts{f};
        if isempty(frontIdx); continue; end
        frontCosts = costMatrix(frontIdx,:);
        cd = computeCrowdingDistance(frontCosts);
        for k = 1:numel(frontIdx)
            pop(frontIdx(k)).Rank = f;
            pop(frontIdx(k)).CrowdingDistance = cd(k);
        end
    end

    % --- Create offspring population ---
    offspring = repmat(empty_ind, nPop, 1);
    for i = 1:2:nPop
        p1 = tournamentSelection(pop);
        p2 = tournamentSelection(pop);

        [c1Pos, c2Pos] = sbxCrossover(p1.Position, p2.Position, ...
            varMin, varMax, etaC, pC);

        c1Pos = polynomialMutation(c1Pos, varMin, varMax, etaM, pM);
        c2Pos = polynomialMutation(c2Pos, varMin, varMax, etaM, pM);

        offspring(i).Position = c1Pos;
        offspring(i).Cost     = objectiveFunction(c1Pos, problem);

        if i+1 <= nPop
            offspring(i+1).Position = c2Pos;
            offspring(i+1).Cost     = objectiveFunction(c2Pos, problem);
        end
    end

    % --- Combine parent and offspring populations ---
    combined       = [pop; offspring];
    costCombined   = vertcat(combined.Cost);
    frontsCombined = nondominatedSort(costCombined);

    % --- Create new population (elitist selection) ---
    newPop  = repmat(empty_ind, nPop, 1);
    nFilled = 0;

    for f = 1:numel(frontsCombined)
        frontIdx = frontsCombined{f};
        if isempty(frontIdx); continue; end
        frontCosts = costCombined(frontIdx,:);
        cd = computeCrowdingDistance(frontCosts);

        if nFilled + numel(frontIdx) <= nPop
            % Take entire front
            for k = 1:numel(frontIdx)
                idx = frontIdx(k);
                nFilled = nFilled + 1;
                newPop(nFilled).Position         = combined(idx).Position;
                newPop(nFilled).Cost             = combined(idx).Cost;
                newPop(nFilled).Rank             = f;
                newPop(nFilled).CrowdingDistance = cd(k);
            end
        else
            % Need only a subset of this front
            [~, sortIdx] = sort(cd, 'descend'); % larger CD first
            nRemain = nPop - nFilled;
            chosen = frontIdx(sortIdx(1:nRemain));
            for k = 1:numel(chosen)
                idx = chosen(k);
                nFilled = nFilled + 1;
                newPop(nFilled).Position         = combined(idx).Position;
                newPop(nFilled).Cost             = combined(idx).Cost;
                newPop(nFilled).Rank             = f;
                newPop(nFilled).CrowdingDistance = cd(sortIdx(k));
            end
            break;
        end
    end

    pop = newPop;

    fprintf('    NSGA-II Generation %d/%d\n', gen, maxGen);
end

% ----- Final Pareto front (only finite / feasible individuals) -----
finalCosts = vertcat(pop.Cost);
finalPos   = vertcat(pop.Position);

finiteIdx  = all(isfinite(finalCosts),2);
finalCosts = finalCosts(finiteIdx,:);
finalPos   = finalPos(finiteIdx,:);

if isempty(finalCosts)
    paretoFront = [];
    paretoSet   = [];
    nsgaRes.population  = pop;
    nsgaRes.paretoFront = paretoFront;
    nsgaRes.paretoSet   = paretoSet;
    return;
end

frontsFinal = nondominatedSort(finalCosts);
F1          = frontsFinal{1};

paretoFront = finalCosts(F1,:);
paretoSet   = finalPos(F1,:);

nsgaRes.population  = pop;
nsgaRes.paretoFront = paretoFront;
nsgaRes.paretoSet   = paretoSet;

end

% -------------------------------------------------------------------------
function ind = tournamentSelection(pop)
%TOURNAMENTSELECTION Binary tournament based on rank + crowding distance.

n = numel(pop);
i1 = randi(n);
i2 = randi(n);
p1 = pop(i1);
p2 = pop(i2);

if p1.Rank < p2.Rank
    ind = p1;
elseif p2.Rank < p1.Rank
    ind = p2;
else
    % Same rank: prefer larger crowding distance
    if p1.CrowdingDistance > p2.CrowdingDistance
        ind = p1;
    else
        ind = p2;
    end
end

end