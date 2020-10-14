function obj = gnccp_obj(X, t, AG, AH)
    N = size(AH,1);
    U = X * ones(N,N) * X';
    F = sum(sum((U.*AG - X * AH * X').^2));
    obj = (1 - abs(t)) * F + t * sum(sum(X.*X));
end