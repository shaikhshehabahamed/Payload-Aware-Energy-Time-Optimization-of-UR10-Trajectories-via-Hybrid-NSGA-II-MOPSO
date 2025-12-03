%% main.m
% Hybrid MOPSO + NSGA-II multi-objective trajectory optimization
% for UR10 with light vs heavy payload scenarios.

clear; clc; close all;
rng(1); % For reproducibility

%% 1. Load base UR10 robot (for limits)
robotBase = loadrobot("universalUR10", ...
    "DataFormat","row", ...
    "Gravity",[0 0 -9.81]);

config  = homeConfiguration(robotBase);
nJoints = numel(config);

% Extract joint limits for moving joints (in order of config)
jointLower = zeros(1,nJoints);
jointUpper = zeros(1,nJoints);
jointIdx   = 0;

for i = 1:numel(robotBase.Bodies)
    jt = robotBase.Bodies{i}.Joint;
    % NOTE: property name is Type (not JointType)
    if ~strcmp(jt.Type,'fixed')
        jointIdx = jointIdx + 1;
        jointLower(jointIdx) = jt.PositionLimits(1);
        jointUpper(jointIdx) = jt.PositionLimits(2);
        if jointIdx == nJoints
            break;
        end
    end
end

%% 2. Define trajectory & constraint settings

nWaypoints  = 5;     % number of joint-space waypoints
nTimeSteps  = 100;   % number of time samples along trajectory

% -------- Relaxed limits so feasible solutions are more likely ----------
velLimit    = 5.0*ones(1,nJoints);    % rad/s  (was 1.5)
torqueLimit = 400*ones(1,nJoints);    % N*m    (was 150)

% Total trajectory time bounds – avoid too aggressive motions
tMin = 5.0;    % seconds (was 2.0)
tMax = 10.0;   % seconds
% -----------------------------------------------------------------------

% Decision variable bounds: [joint waypoints, total time]
varMinAngles = repmat(jointLower, 1, nWaypoints);
varMaxAngles = repmat(jointUpper, 1, nWaypoints);

varMin = [varMinAngles, tMin];
varMax = [varMaxAngles, tMax];

% Base problem template (scenario-independent)
problemTemplate.nJoints      = nJoints;
problemTemplate.nWaypoints   = nWaypoints;
problemTemplate.nVar         = numel(varMin);
problemTemplate.varMin       = varMin;
problemTemplate.varMax       = varMax;
problemTemplate.jointLower   = jointLower;
problemTemplate.jointUpper   = jointUpper;
problemTemplate.velMaxVec    = velLimit;
problemTemplate.torqueMaxVec = torqueLimit;
problemTemplate.tMin         = tMin;
problemTemplate.tMax         = tMax;
problemTemplate.nTimeSteps   = nTimeSteps;

%% 3. Define payload scenarios (light vs heavy)

payloadLight = 1.0;   % kg
payloadHeavy = 10.0;  % kg

% Scenario A: Light payload
robotLight = loadrobot("universalUR10", ...
    "DataFormat","row", ...
    "Gravity",[0 0 -9.81]);
robotLight = attachPayload(robotLight, payloadLight);

problemLight        = problemTemplate;
problemLight.robot  = robotLight;
problemLight.name   = 'Light Payload';

% Scenario B: Heavy payload
robotHeavy = loadrobot("universalUR10", ...
    "DataFormat","row", ...
    "Gravity",[0 0 -9.81]);
robotHeavy = attachPayload(robotHeavy, payloadHeavy);

problemHeavy        = problemTemplate;
problemHeavy.robot  = robotHeavy;
problemHeavy.name   = 'Heavy Payload';

%% 4. Set optimization parameters

% MOPSO parameters (Phase 1)
mopsoParams.nParticles = 40;
mopsoParams.maxIter    = 40;   % G1
mopsoParams.w          = 0.5;
mopsoParams.c1         = 1.5;
mopsoParams.c2         = 2.0;
mopsoParams.archiveMax = 80;

% NSGA-II parameters (Phase 2)
nsga2Params.nPop       = 60;
nsga2Params.maxGen     = 40;   % G2
nsga2Params.pCrossover = 0.9;
nsga2Params.pMutation  = 1/problemTemplate.nVar;
nsga2Params.etaC       = 15;
nsga2Params.etaM       = 20;

%% 5. Run hybrid optimization for both scenarios

fprintf('Running hybrid optimization: Scenario A (Light payload)...\n');
resultsLight = hybridOptimization(problemLight, mopsoParams, nsga2Params);

fprintf('Running hybrid optimization: Scenario B (Heavy payload)...\n');
resultsHeavy = hybridOptimization(problemHeavy, mopsoParams, nsga2Params);

%% 6. Plot Pareto fronts (Energy vs Time) for both scenarios (Hybrid results)

figure; hold on; grid on; box on;

% --- Extract and filter finite points for LIGHT ---
pfLight = resultsLight.Hybrid.paretoFront;
psLight = resultsLight.Hybrid.paretoSet;

finiteLightIdx = all(isfinite(pfLight),2);
pfLight = pfLight(finiteLightIdx,:);
psLight = psLight(finiteLightIdx,:);

% --- Extract and filter finite points for HEAVY ---
pfHeavy = resultsHeavy.Hybrid.paretoFront;
psHeavy = resultsHeavy.Hybrid.paretoSet;

finiteHeavyIdx = all(isfinite(pfHeavy),2);
pfHeavy = pfHeavy(finiteHeavyIdx,:);
psHeavy = psHeavy(finiteHeavyIdx,:);

if ~isempty(pfLight)
    plot(pfLight(:,2), pfLight(:,1), 'o', 'DisplayName','Light - Hybrid');
end
if ~isempty(pfHeavy)
    plot(pfHeavy(:,2), pfHeavy(:,1), 's', 'DisplayName','Heavy - Hybrid');
end

xlabel('Total time [s]');
ylabel('Total energy [arb. units]');
title('Pareto fronts: Hybrid MOPSO+NSGA-II (Light vs Heavy Payload)');
legend('Location','best');

%% 7. Knee point detection for Hybrid solution in each scenario

% Light
if ~isempty(pfLight)
    [kIdxLight, kneePointLight] = detectKneePoint(pfLight);
    kneeDecLight = psLight(kIdxLight,:);
    plot(kneePointLight(2), kneePointLight(1), 'kp', ...
        'MarkerSize',12,'MarkerFaceColor','y', 'DisplayName','Light Knee');
end

% Heavy
if ~isempty(pfHeavy)
    [kIdxHeavy, kneePointHeavy] = detectKneePoint(pfHeavy);
    kneeDecHeavy = psHeavy(kIdxHeavy,:);
    plot(kneePointHeavy(2), kneePointHeavy(1), 'kp', ...
        'MarkerSize',12,'MarkerFaceColor','c', 'DisplayName','Heavy Knee');
end

legend('Location','best');

%% 8. Compute multi-objective metrics table (HV, GD, Spread)

algNames = {'MOPSO'; 'NSGA-II'; 'Hybrid'};

% Scenario A - Light
metricsLight = resultsLight.metrics;
HV_Light     = [metricsLight.MOPSO.HV;     metricsLight.NSGA2.HV;     metricsLight.Hybrid.HV];
GD_Light     = [metricsLight.MOPSO.GD;     metricsLight.NSGA2.GD;     metricsLight.Hybrid.GD];
Spread_Light = [metricsLight.MOPSO.Spread; metricsLight.NSGA2.Spread; metricsLight.Hybrid.Spread];

metricsTableLight = table(algNames, HV_Light, GD_Light, Spread_Light, ...
    'VariableNames', {'Algorithm','HV','GD','Spread'});

disp('Performance metrics - Scenario A (Light payload):');
disp(metricsTableLight);

% Scenario B - Heavy
metricsHeavy = resultsHeavy.metrics;
HV_Heavy     = [metricsHeavy.MOPSO.HV;     metricsHeavy.NSGA2.HV;     metricsHeavy.Hybrid.HV];
GD_Heavy     = [metricsHeavy.MOPSO.GD;     metricsHeavy.NSGA2.GD;     metricsHeavy.Hybrid.GD];
Spread_Heavy = [metricsHeavy.MOPSO.Spread; metricsHeavy.NSGA2.Spread; metricsHeavy.Hybrid.Spread];

metricsTableHeavy = table(algNames, HV_Heavy, GD_Heavy, Spread_Heavy, ...
    'VariableNames', {'Algorithm','HV','GD','Spread'});

disp('Performance metrics - Scenario B (Heavy payload):');
disp(metricsTableHeavy);

%% 9. Plot final optimized knee trajectories

% -------- Scenario A (Light) --------
if exist('kneeDecLight','var') && ~isempty(kneeDecLight)
    [trajLight, isFeasLight] = evaluateRobotTrajectory(kneeDecLight, problemLight);

    if ~isFeasLight
        warning('Knee solution (Light) is infeasible; skipping detailed plots.');
    else
        fprintf('Scenario A (Light) knee solution: Energy = %.4f, Time = %.4f s\n', ...
            trajLight.energyTotal, trajLight.totalTime);

        figure('Name','Knee Trajectory - Light Payload');
        subplot(2,2,1);
        plot(trajLight.t, trajLight.q);
        xlabel('Time [s]'); ylabel('Joint angle [rad]');
        title('Joint Angles (Light)');
        grid on;

        subplot(2,2,2);
        plot(trajLight.t, trajLight.dq);
        xlabel('Time [s]'); ylabel('Joint velocity [rad/s]');
        title('Joint Velocities (Light)');
        grid on;

        subplot(2,2,3);
        plot(trajLight.t, trajLight.tau);
        xlabel('Time [s]'); ylabel('Torque [Nm]');
        title('Joint Torques (Light)');
        grid on;

        subplot(2,2,4);
        plot(trajLight.t, trajLight.powerTotal);
        xlabel('Time [s]'); ylabel('Power [arb.]');
        title(sprintf('Total Power (Light), Energy = %.3f', trajLight.energyTotal));
        grid on;
    end
end

% -------- Scenario B (Heavy) --------
if exist('kneeDecHeavy','var') && ~isempty(kneeDecHeavy)
    [trajHeavy, isFeasHeavy] = evaluateRobotTrajectory(kneeDecHeavy, problemHeavy);

    if ~isFeasHeavy
        warning('Knee solution (Heavy) is infeasible; skipping detailed plots.');
    else
        fprintf('Scenario B (Heavy) knee solution: Energy = %.4f, Time = %.4f s\n', ...
            trajHeavy.energyTotal, trajHeavy.totalTime);

        figure('Name','Knee Trajectory - Heavy Payload');
        subplot(2,2,1);
        plot(trajHeavy.t, trajHeavy.q);
        xlabel('Time [s]'); ylabel('Joint angle [rad]');
        title('Joint Angles (Heavy)');
        grid on;

        subplot(2,2,2);
        plot(trajHeavy.t, trajHeavy.dq);
        xlabel('Time [s]'); ylabel('Joint velocity [rad/s]');
        title('Joint Velocities (Heavy)');
        grid on;

        subplot(2,2,3);
        plot(trajHeavy.t, trajHeavy.tau);
        xlabel('Time [s]'); ylabel('Torque [Nm]');
        title('Joint Torques (Heavy)');
        grid on;

        subplot(2,2,4);
        plot(trajHeavy.t, trajHeavy.powerTotal);
        xlabel('Time [s]'); ylabel('Power [arb.]');
        title(sprintf('Total Power (Heavy), Energy = %.3f', trajHeavy.energyTotal));
        grid on;
    end
end

disp('Done.');