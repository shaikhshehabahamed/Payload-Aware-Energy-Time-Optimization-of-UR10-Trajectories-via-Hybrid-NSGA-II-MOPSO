function mopsoRes = mopsoPhase(problem, params)
%MOPSOPHASE Multi-objective PSO (MOPSO) exploration phase.

nVar   = problem.nVar;
varMin = problem.varMin;
varMax = problem.varMax;

nObj      = 2;  %#ok<NASGU> (for clarity)
nParticles= params.nParticles;
maxIter   = params.maxIter;
w         = params.w;
c1        = params.c1;
c2        = params.c2;
archiveMax= params.archiveMax;

% Particle template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost     = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost     = [];

particle = repmat(empty_particle, nParticles, 1);

% Initialize swarm
for i = 1:nParticles
    particle(i).Position = varMin + rand(1,nVar).*(varMax - varMin);
    particle(i).Velocity = zeros(1,nVar);
    particle(i).Cost     = objectiveFunction(particle(i).Position, problem);
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost     = particle(i).Cost;
end

% Initialize archive
archivePositions = [];
archiveCosts     = [];

% Add initial particles to archive
newPos = vertcat(particle.Position);
newCost= vertcat(particle.Cost);
[archivePositions, archiveCosts] = updateArchive( ...
    archivePositions, archiveCosts, newPos, newCost, archiveMax);

% Main loop
for it = 1:maxIter
    for i = 1:nParticles
        % Select leader from archive (random among archive)
        if ~isempty(archivePositions)
            leaderIdx = randi(size(archivePositions,1));
            gBest     = archivePositions(leaderIdx,:);
        else
            gBest     = particle(i).Best.Position;
        end

        % Update velocity
        r1 = rand(1,nVar);
        r2 = rand(1,nVar);
        particle(i).Velocity = w*particle(i).Velocity ...
            + c1*r1.*(particle(i).Best.Position - particle(i).Position) ...
            + c2*r2.*(gBest - particle(i).Position);

        % Update position
        particle(i).Position = particle(i).Position + particle(i).Velocity;

        % Apply bounds
        particle(i).Position = max(particle(i).Position, varMin);
        particle(i).Position = min(particle(i).Position, varMax);

        % Evaluate
        particle(i).Cost = objectiveFunction(particle(i).Position, problem);

        % Update personal best (multi-objective dominance)
        if dominates(particle(i).Cost, particle(i).Best.Cost)
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost     = particle(i).Cost;
        elseif ~dominates(particle(i).Best.Cost, particle(i).Cost) ...
                && rand < 0.5
            % If neither dominates the other, occasionally update
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost     = particle(i).Cost;
        end
    end

    % Update archive with current swarm
    newPos  = vertcat(particle.Position);
    newCost = vertcat(particle.Cost);
    [archivePositions, archiveCosts] = updateArchive( ...
        archivePositions, archiveCosts, newPos, newCost, archiveMax);

    fprintf('    MOPSO Iteration %d/%d, Archive size = %d\n', ...
        it, maxIter, size(archivePositions,1));
end

% Extract Pareto front from archive
if isempty(archiveCosts)
    paretoFront = [];
    paretoSet   = [];
else
    fronts = nondominatedSort(archiveCosts);
    F1     = fronts{1};
    paretoFront = archiveCosts(F1,:);
    paretoSet   = archivePositions(F1,:);
end

mopsoRes.swarm           = particle;
mopsoRes.archivePositions= archivePositions;
mopsoRes.archiveCosts    = archiveCosts;
mopsoRes.paretoFront     = paretoFront;
mopsoRes.paretoSet       = paretoSet;

end