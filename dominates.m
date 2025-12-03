function flag = dominates(a, b)
%DOMINATES Return true if vector a dominates vector b (minimization).
%   All objectives assumed to be minimized.

flag = all(a <= b) && any(a < b);

end