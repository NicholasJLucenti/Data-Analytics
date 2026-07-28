function [Xi, library_names, L, Y] = run_narmax(X, lambda, lag_order, library_spec, state_names)
%RUN_NARMAX Discover a discrete-time governing map via NARMAX-style sparse
%regression: X_k = Theta([X_{k-1}, ..., X_{k-lag_order}]) * Xi.
%
%   [Xi, library_names, L, Y] = RUN_NARMAX(X, lambda, lag_order, library_spec, state_names)
%
%   Unlike every other variant here, NARMAX never estimates a derivative
%   at all -- not pointwise (standard SINDy), not via integration (weak
%   SINDy). It regresses the CURRENT state directly against a library of
%   candidate functions of the PAST lag_order states. This sidesteps
%   derivative-estimation noise amplification entirely, which is exactly
%   the failure mode that keeps pushing sparse/noisy datasets toward the
%   weak-form route (see selection/select_route.m) -- NARMAX is a
%   genuinely different way of avoiding that problem, not just a
%   noise-robust way of estimating dx/dt.
%
%   IMPLEMENTATION NOTE: the lag_order most recent states are stacked
%   into one "extended state" vector Z = [X_{k-1}, X_{k-2}, ..., X_{k-lag_order}]
%   (D*lag_order columns), and the EXISTING library builder
%   (libraries/build_library.m) is applied to Z exactly as it would be to
%   any ordinary state vector -- this reuses the polynomial and
%   Hill/saturation library machinery unchanged, rather than needing a
%   separate NARMAX-specific library implementation. library_names are
%   then relabeled from the generic 'x<k>' tokens build_library.m
%   produces into readable '<state>_t-<lag>' form.
%
%   Inputs:
%     X            - state matrix (N x D), typically smoothed (NARMAX is
%                    a direct regression on state values, so noise in X
%                    propagates into the fit the same way it would for
%                    standard SINDy -- smoothing still helps)
%     lambda       - sparsity threshold (default 0.1)
%     lag_order    - number of past time steps included in the regressor
%                    (default 1). Larger lag_order captures longer memory
%                    but grows the library size by a factor of lag_order.
%     library_spec - numeric poly_order, or a Hill-flavor struct (see
%                    libraries/build_library.m) -- applied to the STACKED
%                    lag vector, not to X directly
%     state_names  - optional 1xD cell array of channel names, for
%                    readable library_names (default {'x1','x2',...})
%
%   Outputs:
%     Xi            - sparse coefficient matrix (M x D)
%     library_names - cell array of readable term names, e.g.
%                     'mRNA_t-1', 'mRNA_t-1*protein_t-2'
%     L             - the library matrix actually solved against (= Theta(Z)).
%                     Returned so selection/ensemble_sindy.m can bootstrap
%                     over the same system without rebuilding it.
%     Y             - the target matrix actually solved against
%                     (= X(lag_order+1:end, :)). Same rationale as L.

if nargin < 3 || isempty(lag_order)
    lag_order = 1;
end
if nargin < 2 || isempty(lambda)
    lambda = 0.1;
end
if nargin < 4 || isempty(library_spec)
    library_spec = 1;
end

[N, D] = size(X);

if nargin < 5 || isempty(state_names)
    state_names = arrayfun(@(i) sprintf('x%d', i), 1:D, 'UniformOutput', false);
end

if N <= lag_order
    error('run_narmax:tooFewPoints', ...
        'Need more than lag_order (%d) time points (have %d).', lag_order, N);
end

%% 1. Build the stacked lag-history regressor Z and its target Y
numRows = N - lag_order;
Z = zeros(numRows, D * lag_order);
Z_names = cell(1, D * lag_order);

col = 0;
for lag = 1:lag_order
    for d = 1:D
        col = col + 1;
        Z(:, col) = X(lag_order + (1:numRows) - lag, d);
        Z_names{col} = sprintf('%s_t-%d', state_names{d}, lag);
    end
end
Y = X(lag_order+1:end, :);

%% 2. Build the candidate library on Z (reuses existing poly/Hill machinery)
[L, generic_names] = build_library(Z, library_spec);
library_names = local_rename_terms(generic_names, Z_names);

%% 3. Sparse regression
Xi = stlsq_solve(L, Y, lambda, 10);

end


function names_out = local_rename_terms(generic_names, real_names)
    % Replace build_library.m's generic 'x<k>' tokens with readable
    % '<state>_t-<lag>' names. Substitution proceeds from the highest
    % index down so 'x1' can't accidentally clobber a substring inside
    % 'x10', 'x12', etc.
    M = numel(real_names);
    names_out = generic_names;
    for idx = M:-1:1
        pattern = sprintf('x%d(?!\\d)', idx);
        names_out = regexprep(names_out, pattern, real_names{idx});
    end
end