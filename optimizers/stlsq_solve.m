function [Xi, active_mask] = stlsq_solve(Theta, Y, lambda, max_iter)
%STLSQ_SOLVE Generic sequential-thresholded least-squares sparse regression.
%
%   [Xi, active_mask] = STLSQ_SOLVE(Theta, Y, lambda, max_iter)
%
%   Solves Theta * Xi ~= Y by alternating ordinary least squares with
%   hard-thresholding of small coefficients, converging to a sparse
%   coefficient matrix.
%
%   Inputs:
%     Theta    - library/design matrix (K x M). K is "rows of evidence" --
%                for standard SINDy these are time samples; for weak-form
%                SINDy these are test-function integrals. The solver
%                itself doesn't care which.
%     Y        - target matrix (K x D) -- derivatives for standard SINDy,
%                integrated derivative terms for weak-form SINDy.
%     lambda   - hard sparsity threshold (default 0.1)
%     max_iter - number of thresholding iterations (default 10)
%
%   Outputs:
%     Xi          - sparse coefficient matrix (M x D)
%     active_mask - logical (M x D) mask of which terms survived
%
%   This is intentionally variant-agnostic: it's the shared core called by
%   variants/run_standard_sindy.m, variants/run_weak_sindy.m, and
%   selection/ensemble_sindy.m's bootstrap replicates. Keeping it in one
%   place means a change to the sparsification strategy (e.g. swapping in
%   SR3, or a different thresholding schedule) only has to happen once.

if nargin < 4 || isempty(max_iter)
    max_iter = 10;
end
if nargin < 3 || isempty(lambda)
    lambda = 0.1;
end

[K, M] = size(Theta);
if size(Y, 1) ~= K
    error('stlsq_solve:sizeMismatch', ...
        'Theta and Y must have the same number of rows (got %d and %d).', K, size(Y,1));
end
D = size(Y, 2);

Xi = Theta \ Y; % initial unthresholded least-squares guess
active_mask = true(M, D);

for iter = 1:max_iter
    small = abs(Xi) < lambda;
    Xi(small) = 0;
    active_mask = ~small;

    for j = 1:D
        big = active_mask(:, j);
        if any(big)
            Xi(big, j) = Theta(:, big) \ Y(:, j);
        else
            Xi(:, j) = 0; % every candidate term pruned -- no active dynamics found
        end
    end
end

end