function robot = attachPayload(robot, mass)
%ATTACHPAYLOAD Attach a simple payload rigid body to the end-effector.
%   mass: payload mass [kg]

% Simple approximate inertia and CoM
com = [0 0 0.05]; % 5 cm along z from tool
Ixx = 0.001; Iyy = 0.001; Izz = 0.001;
Ixy = 0; Iyz = 0; Ixz = 0;

payload = rigidBody('payload');
jnt = rigidBodyJoint('payloadFixed','fixed');
setFixedTransform(jnt, trvec2tform([0 0 0]));
payload.Joint = jnt;

payload.Mass         = mass;
payload.CenterOfMass = com;
payload.Inertia      = [Ixx Iyy Izz Ixy Iyz Ixz];

eeName = robot.BodyNames{end};
addBody(robot, payload, eeName);

end