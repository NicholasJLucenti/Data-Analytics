function [Xi, library_names, L, Y] = run_standard_sindy(X, dXdt, lambda, poly_order)
%RUN_STANDARD_SINDY Discover governing equations via classic pointwise
%SINDy (finite-difference derivatives + STLSQ).
%
%   [Xi, library_names, L, Y] = RUN_STANDARD_SINDY(X, dXdt, lambda, poly_order)
%
%   Inputs:
%     X          - State matrix (N x D), already cleaned/smoothed upstream
%     dXdt       - Derivative matrix (N x D), already computed upstream
%     lambda     - Sparsity threshold hyperparameter (default 0.1)
%     poly_order - Highest polynomial degree for the candidate library
%                  (default 2)
%
%   Outputs:
%     Xi            - Sparse coefficient matrix (library_size x D)
%     library_names - Cell array describing what each row represents
%     L             - The library matrix actually solved against (= Theta(X)).
%                     Returned so selection/ensemble_sindy.m can bootstrap
%                     over the same system without rebuilding it.
%     Y             - The target matrix actually solved against (= dXdt,
%                     passed straight through). Same rationale as L.
%
%   NOTE: this function does not smooth or differentiate internally --
%   that's the pipeline's job, done once upstream, since weak-form SINDy
%   needs different (in fact, no) derivative preprocessing. Pass in
%   whatever X/dXdt has already been prepared for this route.

if nargin < 4 || isempty(poly_order)
    poly_order = 2;
end
if nargin < 3 || isempty(lambda)
    lambda = 0.1;
end

[N, D] = size(X);
if ~isequal(size(dXdt), [N, D])
    error('run_standard_sindy:sizeMismatch', ...
        'dXdt must be the same size as X (%d x %d).', N, D);
end

[L, library_names] = build_polynomial_library(X, poly_order);
Y = dXdt;

Xi = stlsq_solve(L, Y, lambda, 10);

end