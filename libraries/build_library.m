function [Theta, names] = build_library(X, spec)
%BUILD_LIBRARY Generate a SINDy candidate library -- either a plain
%polynomial library or one of several biologically-motivated
%Hill/saturation flavors -- depending on the type of `spec`.
%
%   [Theta, names] = BUILD_LIBRARY(X, spec)
%
%   If spec is NUMERIC (a scalar poly_order), this is exactly equivalent
%   to calling build_polynomial_library(X, spec) directly. This keeps
%   every existing call site (run_standard_sindy.m, run_weak_sindy.m,
%   run_implicit_sindy.m, benchmarking/simulate_trajectory.m -- all
%   historically called with a bare poly_order) working unchanged.
%
%   If spec is a STRUCT, it must have a .flavor field, one of:
%     'poly_only'             - build_polynomial_library(X, spec.poly_order) only
%     'hill_activation_only'  - constant + linear + same-species Hill activation terms
%     'hill_repression_only'  - constant + linear + same-species Hill repression terms
%     'hill_mixed'            - constant + linear + same-species activation
%                                AND repression + cross-species terms
%     'poly_plus_hill'        - full build_polynomial_library(X, spec.poly_order)
%                                + hill_mixed's same/cross Hill terms
%
%   plus fields (required depending on flavor):
%     .poly_order - used by 'poly_only' and 'poly_plus_hill'
%     .hill_K     - 1xD vector of Hill saturation constants for this
%                   combo (see libraries/estimate_hill_K_candidates.m)
%     .hill_n     - scalar Hill coefficient for this combo
%
%   Output: Theta, names -- exactly like build_polynomial_library.m, just
%   assembled from whichever term sets the flavor calls for.

if isnumeric(spec)
    [Theta, names] = build_polynomial_library(X, spec);
    return
end

if ~isstruct(spec) || ~isfield(spec, 'flavor')
    error('build_library:badSpec', ...
        'spec must be a numeric poly_order, or a struct with a .flavor field.');
end

switch spec.flavor
    case 'poly_only'
        [Theta, names] = build_polynomial_library(X, spec.poly_order);

    case 'hill_activation_only'
        [Theta_lin, names_lin] = build_polynomial_library(X, 1); % constant + linear background
        [Theta_hill, names_hill] = build_hill_library(X, spec.hill_K, spec.hill_n, false);
        keep = startsWith(names_hill, 'hillA');
        Theta = [Theta_lin, Theta_hill(:, keep)];
        names = [names_lin, names_hill(keep)];

    case 'hill_repression_only'
        [Theta_lin, names_lin] = build_polynomial_library(X, 1);
        [Theta_hill, names_hill] = build_hill_library(X, spec.hill_K, spec.hill_n, false);
        keep = startsWith(names_hill, 'hillR');
        Theta = [Theta_lin, Theta_hill(:, keep)];
        names = [names_lin, names_hill(keep)];

    case 'hill_mixed'
        [Theta_lin, names_lin] = build_polynomial_library(X, 1);
        [Theta_hill, names_hill] = build_hill_library(X, spec.hill_K, spec.hill_n, true);
        Theta = [Theta_lin, Theta_hill];
        names = [names_lin, names_hill];

    case 'poly_plus_hill'
        [Theta_poly, names_poly] = build_polynomial_library(X, spec.poly_order);
        [Theta_hill, names_hill] = build_hill_library(X, spec.hill_K, spec.hill_n, true);
        Theta = [Theta_poly, Theta_hill];
        names = [names_poly, names_hill];

    otherwise
        error('build_library:badFlavor', 'Unknown library flavor: %s', spec.flavor);
end

end