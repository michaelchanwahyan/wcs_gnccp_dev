function obj = wcs_obj(X, AG, AH)
    N = size(X,2);
    ONES_NN = ones(N,N);
    U = X * ONES_NN * X';
    obj = sum(sum((U .* AG - X * AH * X').^2));
end