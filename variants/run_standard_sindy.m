function [Xi, library_names] = run_standard_sindy(X, dXdt, lambda, poly_order)
%RUN_STANDARD_SINDY Discover governing equations via classic sequential
%thresholded least squares (STLSQ).
%
%   [Xi, library_names] = RUN_STANDARD_SINDY(X, dXdt, lambda, poly_order)
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
%
%   NOTE: this function no longer smooths or differentiates internally --
%   that is now the shared pipeline's job, done once upstream (see
%   examples/test_hes1_pipeline.m), because weak/implicit SINDy variants
%   need different -- or no -- derivative preprocessing. Pass in whatever
%   X/dXdt has already been prepared for this route.

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

%% 1. Construct candidate library Theta(X)
[Theta, library_names] = build_polynomial_library(X, poly_order);

%% 2. Sequential Thresholded Least Squares (STLSQ)
Xi = Theta \ dXdt; % initial unthresholded least-squares guess

for iter = 1:10
    small_indices = abs(Xi) < lambda;
    Xi(small_indices) = 0;

    for j = 1:D
        big_indices = ~small_indices(:, j);
        if any(big_indices)
            Xi(big_indices, j) = Theta(:, big_indices) \ dXdt(:, j);
        else
            Xi(:, j) = 0; % every term pruned -- no active dynamics found
        end
    end
end

end