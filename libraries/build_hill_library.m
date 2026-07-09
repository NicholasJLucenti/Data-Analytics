function [Theta, names] = build_hill_library(X, K, n, include_cross_species)
%BUILD_HILL_LIBRARY Generate saturating Hill-type candidate terms for
%biological system identification.
%
%   [Theta, names] = BUILD_HILL_LIBRARY(X, K, n, include_cross_species)
%
%   For each state dimension i, adds two same-species terms:
%     activation:  x_i^n / (K_i^n + x_i^n)   -- cooperative binding,
%                  saturating activation, gene induction
%     repression:  K_i^n / (K_i^n + x_i^n)   -- saturating repression,
%                  negative feedback (n=1 is plain Michaelis-Menten form)
%
%   If include_cross_species is true, additionally for every ordered pair
%   i ~= j, adds two cross-species terms:
%     cross-activation: x_i * ( x_j^n / (K_j^n + x_j^n) )
%     cross-repression: x_i * ( K_j^n / (K_j^n + x_j^n) )
%   These represent species i's own (linear/mass-action) term GATED by
%   another species' saturating activation or repression -- e.g. "protein
%   production rate saturating based on mRNA availability", or "mRNA
%   degradation gated by protein-bound repression". This is distinct from
%   (not reducible to) the same-species terms, since it's a genuine
%   product of two different channels.
%
%   Inputs:
%     X - state matrix (N x D)
%     K - 1xD vector of saturation constants, one per channel (typically
%         from libraries/estimate_hill_K_candidates.m)
%     n - scalar Hill coefficient (integer, typically 1, 2, or 4 in
%         biological contexts: non-cooperative, weakly cooperative,
%         strongly cooperative binding respectively)
%     include_cross_species - logical, whether to also add the D*(D-1)*2
%         cross terms above (false gives only the 2*D same-species terms)
%
%   Outputs:
%     Theta - candidate library matrix (N x number_of_terms)
%     names - cell array of human-readable term names, same order as columns

D = size(X, 2);
K = K(:)';
if numel(K) ~= D
    error('build_hill_library:badK', 'K must have one value per channel (numel(K) == size(X,2)).');
end
K(K <= 0) = eps; % guard against a degenerate/zero saturation constant

Theta = [];
names = {};

% Same-species activation/repression
for i = 1:D
    xi_n = X(:, i).^n;
    den_i = K(i)^n + xi_n;
    hillA = xi_n ./ den_i;
    hillR = (K(i)^n) ./ den_i;
    Theta = [Theta, hillA, hillR]; %#ok<AGROW>
    names{end+1} = sprintf('hillA_x%d_n%d_K%.3g', i, n, K(i)); %#ok<AGROW>
    names{end+1} = sprintf('hillR_x%d_n%d_K%.3g', i, n, K(i)); %#ok<AGROW>
end

% Cross-species: species i's own value, gated by species j's saturation
if include_cross_species && D > 1
    for i = 1:D
        for j = 1:D
            if i == j
                continue
            end
            xj_n = X(:, j).^n;
            den_j = K(j)^n + xj_n;
            cross_act = X(:, i) .* (xj_n ./ den_j);
            cross_rep = X(:, i) .* ((K(j)^n) ./ den_j);
            Theta = [Theta, cross_act, cross_rep]; %#ok<AGROW>
            names{end+1} = sprintf('x%d_hillA_x%d_n%d_K%.3g', i, j, n, K(j)); %#ok<AGROW>
            names{end+1} = sprintf('x%d_hillR_x%d_n%d_K%.3g', i, j, n, K(j)); %#ok<AGROW>
        end
    end
end

end