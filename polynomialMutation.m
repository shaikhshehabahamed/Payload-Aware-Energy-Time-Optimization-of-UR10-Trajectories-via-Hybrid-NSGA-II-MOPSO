function y = polynomialMutation(x, varMin, varMax, etaM, pM)
%POLYNOMIALMUTATION Polynomial mutation for real-valued variables.

nVar = numel(x);
y = x;

for i = 1:nVar
    if rand <= pM
        lb = varMin(i);
        ub = varMax(i);
        if ub <= lb
            continue;
        end
        delta1 = (y(i) - lb)/(ub - lb);
        delta2 = (ub - y(i))/(ub - lb);
        rnd = rand;
        mut_pow = 1/(etaM + 1);
        if rnd <= 0.5
            xy = 1 - delta1;
            val = 2*rnd + (1 - 2*rnd)*(xy^(etaM + 1));
            deltaq = val^mut_pow - 1;
        else
            xy = 1 - delta2;
            val = 2*(1 - rnd) + 2*(rnd - 0.5)*(xy^(etaM + 1));
            deltaq = 1 - val^mut_pow;
        end
        y(i) = y(i) + deltaq*(ub - lb);
        % Bound repair
        y(i) = min(max(y(i), lb), ub);
    end
end

end