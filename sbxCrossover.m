function [c1, c2] = sbxCrossover(p1, p2, varMin, varMax, etaC, pC)
%SBXCROSSOVER Simulated Binary Crossover (SBX) for real variables.
%   p1, p2: parent vectors
%   varMin, varMax: bounds (1 x nVar)
%   etaC: distribution index
%   pC: crossover probability

nVar = numel(p1);
c1 = p1;
c2 = p2;

if rand > pC
    return;
end

for i = 1:nVar
    x1 = p1(i);
    x2 = p2(i);
    y1 = min(x1,x2);
    y2 = max(x1,x2);
    lb = varMin(i);
    ub = varMax(i);

    if abs(x1 - x2) > 1e-14
        randu = rand;
        beta = 1 + (2*(y1 - lb)/(y2 - y1));
        alpha = 2 - beta^(-(etaC+1));
        if randu <= 1/alpha
            betaq = (randu*alpha)^(1/(etaC+1));
        else
            betaq = (1/(2 - randu*alpha))^(1/(etaC+1));
        end
        c1_i = 0.5*((y1 + y2) - betaq*(y2 - y1));

        randu = rand;
        beta = 1 + (2*(ub - y2)/(y2 - y1));
        alpha = 2 - beta^(-(etaC+1));
        if randu <= 1/alpha
            betaq = (randu*alpha)^(1/(etaC+1));
        else
            betaq = (1/(2 - randu*alpha))^(1/(etaC+1));
        end
        c2_i = 0.5*((y1 + y2) + betaq*(y2 - y1));

        c1(i) = min(max(c1_i, lb), ub);
        c2(i) = min(max(c2_i, lb), ub);
    else
        c1(i) = x1;
        c2(i) = x2;
    end
end

end