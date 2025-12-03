function f = objectiveFunction(x, problem)
%OBJECTIVEFUNCTION Wrapper around evaluateRobotTrajectory
%   Returns [energy, totalTime]. If any constraint is violated,
%   returns [Inf, Inf].

[traj, isFeasible] = evaluateRobotTrajectory(x, problem);

if ~isFeasible || ~isfinite(traj.energyTotal)
    f = [Inf, Inf];
else
    f = [traj.energyTotal, traj.totalTime];
end

end