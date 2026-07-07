function route = select_route(diagnostics)
%SELECT_ROUTE Decide preprocessing/variant strategy from diagnostics.
%
%   route = SELECT_ROUTE(diagnostics)
%
%   Input:
%     diagnostics - output of diagnostics/profile_data.m (per-channel
%                   struct array plus a .summary rollup)
%
%   Output:
%     route - struct:
%       .recommended_variant  - 'standard' | 'weak'
%       .differentiation      - 'finite_difference' | 'regularized' | 'n/a'
%       .reason               - human-readable explanation
%       .implemented          - true if recommended_variant has a working
%                                variants/*.m script today
%       .warnings             - cell array of caveats (sparse data, low
%                                horizon overlap, etc.)
%
%   This formalizes the triage thresholds that originally lived inline in
%   diagnostics/profile_data.m: SNR < 20 dB -> weak/integral formulation,
%   20-40 dB -> standard SINDy with regularized differentiation, >= 40 dB
%   -> standard SINDy with plain finite differences. Sparsity is now also
%   factored in explicitly: a channel can have a healthy SNR but still be
%   too sparsely sampled (relative to its own dynamical timescale) for
%   safe pointwise derivative estimation.
%
%   Both 'standard' (variants/run_standard_sindy.m) and 'weak'
%   (variants/run_weak_sindy.m) are implemented as of this version.
%   Implicit SINDy is deliberately NOT selectable from diagnostics alone.
%   Whether rational terms are warranted is a structural question about
%   the dynamics, not a noise/sampling question -- it requires fitting a
%   candidate model and testing the residual (see the planned
%   selection/nested_model_test.m), so it is never recommended here.

warnings = {};

minSnr = diagnostics.summary.min_snr_db;
anySparse = diagnostics.summary.any_sparse;
overlapFraction = diagnostics.summary.overlap_fraction;

if overlapFraction < 0.5
    warnings{end+1} = sprintf(...
        'Shared time horizon covers only %.0f%% of the full observed range -- a lot of data is being discarded.', ...
        100*overlapFraction); %#ok<AGROW>
end

if isnan(minSnr)
    warnings{end+1} = 'SNR could not be estimated for at least one channel (too few points); defaulting to the conservative weak-form route.'; %#ok<AGROW>
    route.recommended_variant = 'weak';
    route.differentiation = 'n/a';
    route.reason = 'Insufficient data for noise estimation; weak-form regression is more robust to this uncertainty.';
    route.implemented = true;
    route.warnings = warnings;
    return
end

if anySparse
    warnings{end+1} = 'At least one channel is under-sampled relative to its own dynamical timescale (points_per_period below threshold).'; %#ok<AGROW>
end

if minSnr < 20 || anySparse
    route.recommended_variant = 'weak';
    route.differentiation = 'n/a';
    if minSnr < 20 && anySparse
        route.reason = sprintf(['Worst-case SNR is %.1f dB (< 20 dB) and at least one channel is sparsely ' ...
            'sampled -- pointwise derivative estimation is unreliable. Weak/integral-form SINDy avoids ' ...
            'differentiating the raw signal entirely.'], minSnr);
    elseif minSnr < 20
        route.reason = sprintf(['Worst-case SNR is %.1f dB (< 20 dB). Weak/integral-form SINDy is preferred ' ...
            'because it integrates against test functions instead of differentiating noisy data.'], minSnr);
    else
        route.reason = ['At least one channel is sparsely sampled relative to its dynamics. Weak-form ' ...
            'SINDy is preferred because it does not require a uniform, densely sampled derivative estimate.'];
    end
    route.implemented = true;

elseif minSnr < 40
    route.recommended_variant = 'standard';
    route.differentiation = 'regularized';
    route.reason = sprintf(['Worst-case SNR is %.1f dB (20-40 dB range). Standard SINDy is workable, but ' ...
        'derivatives should ideally use a regularized/robust method rather than plain finite differences.'], minSnr);
    route.implemented = true; % run_standard_sindy.m works; a dedicated regularized-differentiation
                               % scheme is not yet built, so this currently still falls back to
                               % smooth_data + compute_derivatives

else
    route.recommended_variant = 'standard';
    route.differentiation = 'finite_difference';
    route.reason = sprintf(['Worst-case SNR is %.1f dB (>= 40 dB). Data is clean enough for standard ' ...
        'SINDy with plain finite-difference derivatives.'], minSnr);
    route.implemented = true;
end

route.warnings = warnings;

end