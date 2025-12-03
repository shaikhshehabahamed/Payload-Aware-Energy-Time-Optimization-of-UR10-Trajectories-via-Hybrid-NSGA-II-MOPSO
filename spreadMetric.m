function delta = spreadMetric(front, refFront)
%SPREADMETRIC Deb's spread (Δ) metric.
%   front: N x 2
%   refFront: K x 2 (reference front)

if size(front,1) < 2 || size(refFront,1) < 2
    delta = Inf;
    return;
end

% Sort both fronts by first objective
frontSorted = sortrows(front,1);
refSorted   = sortrows(refFront,1);

e1 = refSorted(1,:);
e2 = refSorted(end,:);

df = norm(frontSorted(1,:)   - e1);
dl = norm(frontSorted(end,:) - e2);

% Distances between consecutive points
N = size(frontSorted,1);
di = zeros(N-1,1);
for i = 1:(N-1)
    di(i) = norm(frontSorted(i+1,:) - frontSorted(i,:));
end

dbar = mean(di);
if dbar == 0
    delta = Inf;
    return;
end

delta = (df + dl + sum(abs(di - dbar))) / (df + dl + (N-1)*dbar + eps);

end