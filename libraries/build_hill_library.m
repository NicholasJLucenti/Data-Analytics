function [Theta, names] = build_hill_library(X, K, n, include_cross_species, same_forms, cross_forms)
%BUILD_HILL_LIBRARY Generate saturating Hill-type candidate terms for
%biological system identification.
%
%   [Theta, names] = BUILD_HILL_LIBRARY(X, K, n, include_cross_species, same_forms, cross_forms)
%
%   For each state dimension i, can add same-species terms:
%     activation:  x_i^n / (K_i^n + x_i^n)   -- cooperative binding,
%                  saturating activation, gene induction
%     repression:  K_i^n / (K_i^n + x_i^n)   -- saturating repression,
%                  negative feedback (n=1 is plain Michaelis-Menten form)
%
%   And, if include_cross_species, for every ordered pair i ~= j:
%     cross-activation: x_i * ( x_j^n / (K_j^n + x_j^n) )
%     cross-repression: x_i * ( K_j^n / (K_j^n + x_j^n) )
%
%   IMPORTANT -- exact linear redundancy: activation_i + repression_i = 1
%   EXACTLY (they're complementary fractions of the same denominator), and
%   cross_act_i_j + cross_rep_i_j = x_i EXACTLY. Any library that already
%   contains a constant term (activation+repression case) or the linear
%   term x_i (cross case) -- which every flavor in build_library.m does --
%   becomes rank-deficient if BOTH forms of the SAME target are included
%   together, since one is then an exact linear combination of the other
%   two. same_forms/cross_forms let the caller select only ONE form per
%   target to avoid this. Requesting both forms for the same target
%   (e.g. same_forms={'activation','repression'}) is only safe when no
%   constant/linear background term is present elsewhere in the library
%   -- build_library.m never does this for you automatically, so pick
%   deliberately.
%
%   Inputs:
%     X    - state matrix (N x D)
%     K    - 1xD vector of saturation constants, one per channel
%     n    - scalar Hill coefficient
%     include_cross_species - logical, whether to add cross-species terms
%     same_forms  - cell array subset of {'activation','repression'},
%                   which same-species forms to generate (default: both --
%                   only safe if you know what you're doing, see above)
%     cross_forms - cell array subset of {'activation','repression'},
%                   which cross-species forms to generate (default: both,
%                   same caveat; ignored if include_cross_species is false)
%
%   Outputs:
%     Theta - candidate library matrix (N x number_of_terms)
%     names - cell array of human-readable term names, same order as columns

if nargin < 4 || isempty(include_cross_species)
    include_cross_species = true;
end
if nargin < 5 || isempty(same_forms)
    same_forms = {'activation', 'repression'};
end
if nargin < 6 || isempty(cross_forms)
    cross_forms = {'activation', 'repression'};
end

D = size(X, 2);
K = K(:)';
if numel(K) ~= D
    error('build_hill_library:badK', 'K must have one value per channel (numel(K) == size(X,2)).');
end
K(K <= 0) = eps; % guard against a degenerate/zero saturation constant

Theta = [];
names = {};

% Same-species terms
for i = 1:D
    xi_n = X(:, i).^n;
    den_i = K(i)^n + xi_n;
    if ismember('activation', same_forms)
        Theta = [Theta, xi_n ./ den_i]; %#ok<AGROW>
        names{end+1} = sprintf('hillA_x%d_n%d_K%.3g', i, n, K(i)); %#ok<AGROW>
    end
    if ismember('repression', same_forms)
        Theta = [Theta, (K(i)^n) ./ den_i]; %#ok<AGROW>
        names{end+1} = sprintf('hillR_x%d_n%d_K%.3g', i, n, K(i)); %#ok<AGROW>
    end
end

% Cross-species terms: species i's own value, gated by species j's saturation
if include_cross_species && D > 1
    for i = 1:D
        for j = 1:D
            if i == j
                continue
            end
            xj_n = X(:, j).^n;
            den_j = K(j)^n + xj_n;
            if ismember('activation', cross_forms)
                Theta = [Theta, X(:, i) .* (xj_n ./ den_j)]; %#ok<AGROW>
                names{end+1} = sprintf('x%d_hillA_x%d_n%d_K%.3g', i, j, n, K(j)); %#ok<AGROW>
            end
            if ismember('repression', cross_forms)
                Theta = [Theta, X(:, i) .* ((K(j)^n) ./ den_j)]; %#ok<AGROW>
                names{end+1} = sprintf('x%d_hillR_x%d_n%d_K%.3g', i, j, n, K(j)); %#ok<AGROW>
            end
        end
    end
end

end