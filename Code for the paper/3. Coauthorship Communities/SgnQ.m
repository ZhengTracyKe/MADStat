function [Q_pv,Qn,Q_st] = SgnQ(A)
A = A - diag(diag(A));
[n,~] = size(A);
V = ones(1,n) * A * ones(n,1);
eta = 1 / V^0.5 * A * ones(n,1);

At = A - eta * eta';

dAt = diag(diag(At));

Qn = trace(At^4) - 4 * trace(At .* At^3) + 8 * trace(At .* At .* At^2) ...
    - 6 * trace(At .* At .* At .* At) - 2 * trace(At^2 .* At^2) ...
    + 2 * ones(1,n) * dAt * (At .* At) * dAt * ones(n,1) ...
    + ones(1,n) * (At .* At .* At .* At) * ones(n,1);

Q_st = (Qn - 2 * (sum(eta.^2) - 1)^2) / (8^0.5 * (sum(eta.^2) - 1)^2);

Q_pv = 1 - normcdf(Q_st);

end