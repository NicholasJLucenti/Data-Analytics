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
%     Theta    - library/design matrix (K x M)
%     Y        - target matrix (K x D)
%     lambda   - hard sparsity threshold (default 0.1)
%     max_iter - number of thresholding iterations (default 10)
%
%   Outputs:
%     Xi          - sparse coefficient matrix (M x D), in Theta's ORIGINAL units
%     active_mask - logical (M x D) mask of which terms survived
%
%   COLUMN NORMALIZATION: lambda is a single global cutoff, so without
%   normalization it implicitly favors/penalizes terms based on their raw
%   numeric scale rather than their actual contribution. This matters a
%   lot for run_implicit_sindy.m, whose design matrix mixes numerator
%   terms (~x^k) with interaction terms (~x^k * dxdt) that can differ by
%   orders of magnitude -- a uniform threshold on unnormalized
%   coefficients can leave a poorly-conditioned (near-zero-crossing)
%   denominator purely as a scale artifact. Each column of Theta is
%   normalized to unit L2 norm before thresholding, and coefficients are
%   rescaled back to Theta's original units once at the end -- so lambda
%   always means "this term's normalized-scale contribution is small",
%   consistently across every variant that calls this solver.
%
%   Uses lsqminnorm rather than backslash (\) for every least-squares
%   solve, so a rank-deficient or ill-conditioned Theta degrades
%   gracefully (well-defined minimum-norm solution) instead of printing a
%   warning and returning an unpredictable result.
%
%   This is intentionally variant-agnostic: it's the shared core called by
%   variants/run_standard_sindy.m, variants/run_weak_sindy.m,
%   variants/run_implicit_sindy.m, and selection/ensemble_sindy.m's
%   bootstrap replicates.

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

colNorms = sqrt(sum(Theta.^2, 1));
colNorms(colNorms == 0) = 1; % guard an all-zero column (shouldn't occur in practice)
Theta_n = Theta ./ colNorms;

Xi_n = lsqminnorm(Theta_n, Y); % initial unthresholded least-squares guess, normalized scale
active_mask = true(M, D);

for iter = 1:max_iter
    small = abs(Xi_n) < lambda;
    Xi_n(small) = 0;
    active_mask = ~small;

    for j = 1:D
        big = active_mask(:, j);
        if any(big)
            Xi_n(big, j) = lsqminnorm(Theta_n(:, big), Y(:, j));
        else
            Xi_n(:, j) = 0; % every candidate term pruned -- no active dynamics found
        end
    end
end

Xi = Xi_n ./ colNorms'; % rescale back to Theta's original units

end