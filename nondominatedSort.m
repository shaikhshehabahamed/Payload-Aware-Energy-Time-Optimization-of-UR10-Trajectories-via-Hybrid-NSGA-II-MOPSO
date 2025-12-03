function fronts = nondominatedSort(costs)
%NONDOMINATEDSORT Fast non-dominated sorting (Deb's NSGA-II method)
%   costs: N x M matrix of objective values (minimization)
%   fronts: cell array, each element is vector of indices

N = size(costs,1);
S = cell(N,1);
nDom = zeros(N,1);
fronts = {};

F1 = [];

for p = 1:N
    S{p} = [];
    nDom(p) = 0;
    for q = 1:N
        if p == q
            continue;
        end
        if dominates(costs(p,:), costs(q,:))
            S{p} = [S{p}, q]; %#ok<AGROW>
        elseif dominates(costs(q,:), costs(p,:))
            nDom(p) = nDom(p) + 1;
        end
    end
    if nDom(p) == 0
        F1 = [F1, p]; %#ok<AGROW>
    end
end

fronts{1} = F1;
i = 1;

while ~isempty(fronts{i})
    Q = [];
    for p = fronts{i}
        for q = S{p}
            nDom(q) = nDom(q) - 1;
            if nDom(q) == 0
                Q = [Q, q]; %#ok<AGROW>
            end
        end
    end
    i = i + 1;
    fronts{i} = unique(Q);
end

% Remove possible last empty front
if isempty(fronts{end})
    fronts(end) = [];
end

end