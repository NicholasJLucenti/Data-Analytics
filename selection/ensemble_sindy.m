function ensemble = ensemble_sindy(L, Y, lambda, varargin)
%ENSEMBLE_SINDY Bootstrap-aggregate (bagging) sparse regression over a
%SINDy library/target system to get coefficient uncertainty instead of a
%single point estimate.
%
%   ensemble = ENSEMBLE_SINDY(L, Y, lambda, ...)
%
%   Rather than solving Theta*Xi = Y (or the weak-form A*xi = b) once,
%   this resamples rows of the system with replacement many times, refits
%   STLSQ on each resample, and aggregates. Terms that are real dynamics
%   should be selected consistently across resamples; terms that are
%   noise artifacts of the particular dataset will drop in and out. This
%   directly answers "how much do I trust this coefficient" -- which a
%   single STLSQ point estimate cannot.
%
%   This function is variant-agnostic: it operates purely on whatever
%   (L, Y) linear system a variant produces, not on raw time-series data.
%   That means it works unchanged for:
%     - standard SINDy: L = Theta(X), Y = dXdt  (rows = time samples)
%     - weak SINDy:     L = A,        Y = b     (rows = test functions)
%   Get L and Y from either variant's 3rd/4th output
%   (e.g. [~, ~, L, Y] = run_standard_sindy(...) or run_weak_sindy(...)),
%   then pass them here.
%
%   Inputs:
%     L      - library/design matrix (K x M)
%     Y      - target matrix (K x D)
%     lambda - sparsity threshold used in every bootstrap refit
%
%   Name-value options:
%     'NumReplicates'     - number of bootstrap resamples (default 200)
%     'SampleFraction'    - resample size as a fraction of K, sampled
%                           WITH replacement (default 1.0, i.e. a
%                           classic same-size bootstrap sample)
%     'InclusionThreshold' - a term is kept in the final aggregated model
%                           only if it was nonzero in at least this
%                           fraction of replicates (default 0.5, i.e.
%                           majority vote). This inclusion-probability
%                           mechanism is also the natural hook for the
%                           planned HPO grid search over sparsity/
%                           inclusion parameters -- inclusion_prob is
%                           exactly the quantity that grid search would
%                           be tuning against.
%     'Seed'              - RNG seed for reproducibility (default: not set)
%
%   Output: ensemble struct with fields:
%     .Xi_final       - (M x D) aggregated sparse model: for each term
%                       that clears InclusionThreshold, the median of its
%                       nonzero bootstrap estimates; zero elsewhere.
%     .Xi_mean        - (M x D) mean coefficient across all replicates
%                       (zeros included, so this shrinks toward zero for
%                       inconsistently-selected terms)
%     .Xi_median      - (M x D) median coefficient across all replicates
%     .Xi_std         - (M x D) std of coefficient across all replicates
%     .inclusion_prob - (M x D) fraction of replicates in which each term
%                       was nonzero -- the key stability diagnostic
%     .Xi_stack       - (NumReplicates x M x D) raw per-replicate results,
%                       kept for any further custom analysis
%     .num_replicates, .sample_size, .lambda - bookkeeping

p = inputParser;
addParameter(p, 'NumReplicates', 200, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'SampleFraction', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'InclusionThreshold', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(p, 'Seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
opts = p.Results;

if nargin < 3 || isempty(lambda)
    lambda = 0.1;
end

[K, M] = size(L);
if size(Y, 1) ~= K
    error('ensemble_sindy:sizeMismatch', 'L and Y must have the same number of rows.');
end
D = size(Y, 2);

if ~isempty(opts.Seed)
    rng(opts.Seed);
end

numReplicates = round(opts.NumReplicates);
sampleSize = max(round(opts.SampleFraction * K), 1);

Xi_stack = zeros(numReplicates, M, D);

for r = 1:numReplicates
    idx = randi(K, sampleSize, 1); % bootstrap resample with replacement
    Lr = L(idx, :);
    Yr = Y(idx, :);
    Xi_stack(r, :, :) = stlsq_solve(Lr, Yr, lambda, 10);
end

Xi_mean = reshape(mean(Xi_stack, 1), [M, D]);
Xi_median = reshape(median(Xi_stack, 1), [M, D]);
Xi_std = reshape(std(Xi_stack, 0, 1), [M, D]);
inclusion_prob = reshape(mean(Xi_stack ~= 0, 1), [M, D]);

Xi_final = zeros(M, D);
for d = 1:D
    for m = 1:M
        if inclusion_prob(m, d) >= opts.InclusionThreshold
            vals = Xi_stack(:, m, d);
            nz = vals(vals ~= 0);
            if ~isempty(nz)
                Xi_final(m, d) = median(nz);
            end
        end
    end
end

ensemble.Xi_final = Xi_final;
ensemble.Xi_mean = Xi_mean;
ensemble.Xi_median = Xi_median;
ensemble.Xi_std = Xi_std;
ensemble.inclusion_prob = inclusion_prob;
ensemble.Xi_stack = Xi_stack;
ensemble.num_replicates = numReplicates;
ensemble.sample_size = sampleSize;
ensemble.lambda = lambda;

end