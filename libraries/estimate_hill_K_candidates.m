function K_candidates = estimate_hill_K_candidates(X, num_candidates)
%ESTIMATE_HILL_K_CANDIDATES Derive a small, data-grounded set of
%candidate Hill/Michaelis-Menten saturation constants K, one full K
%vector (one value per channel) per candidate row.
%
%   K_candidates = ESTIMATE_HILL_K_CANDIDATES(X, num_candidates)
%
%   Rather than treating K as a free continuous hyperparameter (which
%   would blow up the grid search combinatorics), K is derived from
%   quantiles of each channel's own observed range -- a biologically
%   meaningful saturation constant is, almost by construction, some
%   fraction of the variable's typical operating range.
%
%   Inputs:
%     X              - state matrix (N x D)
%     num_candidates - number of K choices to generate (default 3).
%                      Kept deliberately small -- this directly controls
%                      how many extra grid-search combinations the Hill
%                      library flavors add.
%
%   Output:
%     K_candidates - (num_candidates x D) matrix. K_candidates(k, :) is
%                    one complete K vector to try as a single combo.
%                    Default (num_candidates=3): 25th/50th/75th
%                    percentiles of each channel.

if nargin < 2 || isempty(num_candidates)
    num_candidates = 3;
end

D = size(X, 2);

if num_candidates == 1
    quantile_levels = 0.5;
else
    quantile_levels = linspace(0.25, 0.75, num_candidates);
end

K_candidates = zeros(num_candidates, D);
for d = 1:D
    K_candidates(:, d) = quantile(X(:, d), quantile_levels)';
end

K_candidates(K_candidates <= 0) = eps; % guard a degenerate/flat channel

end