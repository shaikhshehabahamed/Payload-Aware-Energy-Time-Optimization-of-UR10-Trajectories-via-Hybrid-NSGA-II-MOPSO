function [traj, isFeasible] = evaluateRobotTrajectory(x, problem)
%EVALUATEROBOTTRAJECTORY Build joint spline trajectory, compute dynamics
%   x: decision vector [q_waypoints (nWaypoints*nJoints), totalTime]

nJoints    = problem.nJoints;
nWaypoints = problem.nWaypoints;
nTimeSteps = problem.nTimeSteps;

jointLower = problem.jointLower;
jointUpper = problem.jointUpper;
velMax     = problem.velMaxVec;
torqueMax  = problem.torqueMaxVec;

robot = problem.robot;

% Decode decision vector
T = x(end);
traj.totalTime = T;

% Time bounds check
if T < problem.tMin || T > problem.tMax
    isFeasible = false;
    traj.energyTotal = Inf;
    traj.t = [];
    traj.q = [];
    traj.dq = [];
    traj.ddq = [];
    traj.tau = [];
    traj.powerTotal = [];
    return;
end

qWaypoints = zeros(nWaypoints, nJoints);
for k = 1:nWaypoints
    idxStart = (k-1)*nJoints + 1;
    idxEnd   = k*nJoints;
    qWaypoints(k,:) = x(idxStart:idxEnd);
end
traj.qWaypoints = qWaypoints;

% Joint limit check at waypoints
for j = 1:nJoints
    if any(qWaypoints(:,j) < jointLower(j) - 1e-8) || ...
       any(qWaypoints(:,j) > jointUpper(j) + 1e-8)
        isFeasible = false;
        traj.energyTotal = Inf;
        traj.t = [];
        traj.q = [];
        traj.dq = [];
        traj.ddq = [];
        traj.tau = [];
        traj.powerTotal = [];
        return;
    end
end

% Time discretization
t = linspace(0, T, nTimeSteps);
wpTimes = linspace(0, T, nWaypoints);
traj.t      = t;
traj.wpTimes= wpTimes;

% Build B-spline (cubic spline) joint trajectories
q   = zeros(nTimeSteps, nJoints);
dq  = zeros(nTimeSteps, nJoints);
ddq = zeros(nTimeSteps, nJoints);

for j = 1:nJoints
    pp   = spline(wpTimes, qWaypoints(:,j)');
    q(:,j)   = ppval(pp, t);
    pp_d1 = fnder(pp,1);
    dq(:,j) = ppval(pp_d1, t);
    pp_d2 = fnder(pp,2);
    ddq(:,j)= ppval(pp_d2, t);
end

% Joint limit check over full trajectory
for j = 1:nJoints
    if any(q(:,j) < jointLower(j) - 1e-6) || ...
       any(q(:,j) > jointUpper(j) + 1e-6)
        isFeasible = false;
        traj.energyTotal = Inf;
        traj.q = q; traj.dq = dq; traj.ddq = ddq;
        traj.tau = []; traj.powerTotal = [];
        return;
    end
end

% Velocity limit check
for j = 1:nJoints
    if any(abs(dq(:,j)) > velMax(j) + 1e-8)
        isFeasible = false;
        traj.energyTotal = Inf;
        traj.q = q; traj.dq = dq; traj.ddq = ddq;
        traj.tau = []; traj.powerTotal = [];
        return;
    end
end

% Dynamics: inverseDynamics to get torque, then power and energy
tau = zeros(nTimeSteps, nJoints);
for i = 1:nTimeSteps
    tau(i,:) = inverseDynamics(robot, q(i,:), dq(i,:), ddq(i,:));
end

% Torque limit check
for j = 1:nJoints
    if any(abs(tau(:,j)) > torqueMax(j) + 1e-8)
        isFeasible = false;
        traj.energyTotal = Inf;
        traj.q = q; traj.dq = dq; traj.ddq = ddq;
        traj.tau = tau; traj.powerTotal = [];
        return;
    end
end

% Instantaneous power (per joint) and total absolute power
powerJoint = tau .* dq;
powerTotal = sum(abs(powerJoint), 2);

% Integrate to get total energy (approx)
energyTotal = trapz(t, powerTotal);

isFeasible = true;
traj.q     = q;
traj.dq    = dq;
traj.ddq   = ddq;
traj.tau   = tau;
traj.powerJoint = powerJoint;
traj.powerTotal = powerTotal;
traj.energyTotal= energyTotal;

end