function [results, t, X, names] = run_grid_search(raw, variant, varargin)
%RUN_GRID_SEARCH Sweep SINDy hyperparameters -- including candidate
%library flavor -- and score every resulting model, so a downstream
%selection step (select_best_model.m) can pick among competing
%trade-offs instead of committing to one fixed-parameter fit.
%
%   [results, t, X, names] = RUN_GRID_SEARCH(raw, variant, ...)
%
%   Inputs:
%     raw     - raw struct array (io/load_raw_data.m format)
%     variant - 'standard', 'weak', 'implicit', or 'narmax'
%
%   Name-value options:
%     'Resolution'            - 'fast' | 'balanced' (default) | 'thorough'.
%                                Sets ALL the grids below at once via
%                                get_resolution_preset.m. Any grid you
%                                also pass explicitly overrides just that
%                                one field on top of the preset.
%     'LambdaGrid'            - sparsity thresholds to sweep
%     'PolyOrderGrid'         - candidate library orders to sweep
%     'WindowPointsGrid'      - (weak only) test-function window sizes
%     'TestFunctionOrderGrid' - (weak only) test-function smoothness orders
%     'LagOrderGrid'          - (narmax only) number of past time steps
%                                included in the regressor
%     'NumDensePoints'        - densification target
%     'SmoothingFactor'       - (standard/implicit/narmax only) Lowess span
%     'LibraryFlavorGrid'     - cell array of library flavors to sweep,
%                                for standard/weak variants only (implicit
%                                and narmax always use 'poly_only' -- see
%                                notes below). Available: 'poly_only',
%                                'hill_activation_only',
%                                'hill_repression_only', 'hill_mixed',
%                                'poly_plus_hill'.
%     'HillCoeffGrid'         - Hill exponents n to sweep for any hill_*
%                                flavor
%     'HillKCandidates'       - number of data-derived K candidate rows
%
%   See get_resolution_preset.m for what each Resolution level sets, and
%   approximate combination counts per level.
%
%   IMPLICIT VARIANT NOTE: run_implicit_sindy.m's own algebra already
%   produces rational/saturating output from a plain polynomial library;
%   layering Hill terms underneath an already-implicit formulation adds
%   interpretive complexity without a clear benefit, so the implicit
%   variant always sweeps 'poly_only' regardless of LibraryFlavorGrid.
%
%   NARMAX VARIANT NOTE: NARMAX applies build_library.m to a STACKED
%   lag-history vector (see run_narmax.m), so a Hill flavor there would
%   mean "saturating function of a past-lagged value" rather than "of the
%   current state" -- a valid but unusual construction that also
%   multiplies combinatorics by D*lag_order instead of D. NARMAX always
%   sweeps 'poly_only'; LagOrderGrid is its distinguishing sweep
%   dimension instead.
%
%   Output:
%     results - struct array, one entry per parameter combination:
%                 .variant, .flavor
%                 .lambda, .poly_order (NaN for hill_activation_only/
%                   hill_repression_only/hill_mixed, which don't use it),
%                   .window_points, .test_function_order (NaN unless variant=='weak'),
%                   .lag_order (NaN unless variant=='narmax')
%                 .hill_n, .hill_K (NaN/[] for 'poly_only')
%                 .library_spec      - the EXACT spec (numeric poly_order,
%                   or Hill-flavor struct) passed to the variant function.
%                   Use this (not .poly_order) to reconstruct/re-simulate
%                   the model -- the other fields above are for display only.
%                 .num_active_terms, .trajectory_rmse, .normalized_rmse,
%                   .dynamics_class, .period_error, .Xi, .library_names, .success
%     t, X, names - the shared densified trajectory/channel names

[resolution, varargin] = local_extract_option(varargin, 'Resolution', 'balanced');
preset = get_resolution_preset(resolution);

p = inputParser;
addParameter(p, 'LambdaGrid', preset.LambdaGrid, @isnumeric);
addParameter(p, 'PolyOrderGrid', preset.PolyOrderGrid, @isnumeric);
addParameter(p, 'WindowPointsGrid', preset.WindowPointsGrid, @isnumeric);
addParameter(p, 'TestFunctionOrderGrid', preset.TestFunctionOrderGrid, @isnumeric);
addParameter(p, 'LagOrderGrid', preset.LagOrderGrid, @isnumeric);
addParameter(p, 'NumDensePoints', preset.NumDensePoints, @isnumeric);
addParameter(p, 'SmoothingFactor', 0.2, @isnumeric);
addParameter(p, 'LibraryFlavorGrid', preset.LibraryFlavorGrid, @iscell);
addParameter(p, 'HillCoeffGrid', preset.HillCoeffGrid, @isnumeric);
addParameter(p, 'HillKCandidates', preset.HillKCandidates, @isnumeric);
parse(p, varargin{:});
opts = p.Results;

if ~ismember(variant, {'standard', 'weak', 'implicit', 'narmax'})
    error('run_grid_search:badVariant', 'variant must be ''standard'', ''weak'', ''implicit'', or ''narmax''.');
end

%% Shared preprocessing -- done once, independent of the swept parameters
[t, X, names] = align_and_truncate(raw, opts.NumDensePoints);

needsSmoothing = strcmp(variant, 'standard') || strcmp(variant, 'implicit') || strcmp(variant, 'narmax');
needsDerivatives = strcmp(variant, 'standard') || strcmp(variant, 'implicit');
if needsSmoothing
    X_smooth = smooth_data(t, X, opts.SmoothingFactor);
end
if needsDerivatives
    dXdt = compute_derivatives(t, X_smooth);
end

K_candidates = estimate_hill_K_candidates(X, opts.HillKCandidates); % (numK x D), data-derived
numK = size(K_candidates, 1);

real_period = estimate_trajectory_period(t, X); % shared target for every combo's period_error below

if strcmp(variant, 'implicit') || strcmp(variant, 'narmax')
    flavorsToSweep = {'poly_only'};
else
    flavorsToSweep = opts.LibraryFlavorGrid;
end

if strcmp(variant, 'weak')
    windowList = opts.WindowPointsGrid;
    tfoList = opts.TestFunctionOrderGrid;
else
    windowList = NaN;
    tfoList = NaN;
end

if strcmp(variant, 'narmax')
    lagList = opts.LagOrderGrid;
else
    lagList = NaN;
end

%% Build the flat list of runs -- different flavors/variants sweep
%% different sub-grids, so this is assembled explicitly per flavor rather
%% than forced into one rectangular ndgrid.
runs = struct('lambda', {}, 'library_spec', {}, 'window_points', {}, 'test_function_order', {}, ...
    'lag_order', {}, 'flavor', {});

for f = 1:numel(flavorsToSweep)
    flavor = flavorsToSweep{f};

    switch flavor
        case 'poly_only'
            for lam = opts.LambdaGrid
                for po = opts.PolyOrderGrid
                    for wp = windowList
                        for tf = tfoList
                            for lg = lagList
                                runs(end+1) = struct('lambda', lam, 'library_spec', po, ...
                                    'window_points', wp, 'test_function_order', tf, ...
                                    'lag_order', lg, 'flavor', flavor); %#ok<AGROW>
                            end
                        end
                    end
                end
            end

        case {'hill_activation_only', 'hill_repression_only', 'hill_mixed'}
            for lam = opts.LambdaGrid
                for hn = opts.HillCoeffGrid
                    for ki = 1:numK
                        spec = struct('flavor', flavor, 'hill_K', K_candidates(ki, :), 'hill_n', hn);
                        for wp = windowList
                            for tf = tfoList
                                runs(end+1) = struct('lambda', lam, 'library_spec', spec, ...
                                    'window_points', wp, 'test_function_order', tf, ...
                                    'lag_order', NaN, 'flavor', flavor); %#ok<AGROW>
                            end
                        end
                    end
                end
            end

        case 'poly_plus_hill'
            for lam = opts.LambdaGrid
                for po = opts.PolyOrderGrid
                    for hn = opts.HillCoeffGrid
                        for ki = 1:numK
                            spec = struct('flavor', flavor, 'poly_order', po, ...
                                'hill_K', K_candidates(ki, :), 'hill_n', hn);
                            for wp = windowList
                                for tf = tfoList
                                    runs(end+1) = struct('lambda', lam, 'library_spec', spec, ...
                                        'window_points', wp, 'test_function_order', tf, ...
                                        'lag_order', NaN, 'flavor', flavor); %#ok<AGROW>
                                end
                            end
                        end
                    end
                end
            end

        otherwise
            error('run_grid_search:badFlavor', 'Unknown flavor in LibraryFlavorGrid: %s', flavor);
    end
end

nCombos = numel(runs);
fprintf('[GRID SEARCH] Sweeping %d parameter combinations (%s variant, resolution=%s, flavors: %s)...\n', ...
    nCombos, variant, resolution, strjoin(flavorsToSweep, ', '));

results = struct('variant', {}, 'flavor', {}, 'lambda', {}, 'poly_order', {}, 'window_points', {}, ...
    'test_function_order', {}, 'lag_order', {}, 'hill_n', {}, 'hill_K', {}, 'library_spec', {}, ...
    'num_active_terms', {}, 'trajectory_rmse', {}, 'normalized_rmse', {}, 'dynamics_class', {}, ...
    'period_error', {}, 'Xi', {}, 'library_names', {}, 'success', {});

progressStep = max(1, round(nCombos / 10));

for c = 1:nCombos
    run_c = runs(c);
    lambda = run_c.lambda;
    spec = run_c.library_spec; % numeric poly_order, OR a Hill-flavor struct
    window_points = run_c.window_points;
    test_function_order = run_c.test_function_order;
    lag_order = run_c.lag_order;

    if isnumeric(spec)
        poly_order_field = spec;
        hill_n_field = NaN;
        hill_K_field = [];
    else
        poly_order_field = local_getfield_default(spec, 'poly_order', NaN);
        hill_n_field = spec.hill_n;
        hill_K_field = spec.hill_K;
    end

    try
        switch variant
            case 'standard'
                [model, library_names] = run_standard_sindy(X_smooth, dXdt, lambda, spec);
                num_active = nnz(model);
            case 'weak'
                [model, library_names] = run_weak_sindy(t, X, lambda, spec, ...
                    'WindowPoints', window_points, 'TestFunctionOrder', test_function_order);
                num_active = nnz(model);
            case 'implicit'
                [model, library_names] = run_implicit_sindy(X_smooth, dXdt, lambda, spec);
                num_active = 0;
                for d = 1:numel(model)
                    num_active = num_active + nnz(model(d).numerator_Xi) + nnz(model(d).denominator_Xi(2:end));
                end
            case 'narmax'
                [model, library_names] = run_narmax(X_smooth, lambda, lag_order, spec, names);
                num_active = nnz(model);
        end

        if strcmp(variant, 'narmax')
            metrics = compute_narmax_trajectory_error(t, X, model, lag_order, spec);
        else
            metrics = compute_trajectory_error(t, X, model, spec);
        end

        if metrics.success
            if strcmp(variant, 'narmax')
                [t_sim, X_sim, sim_info] = simulate_narmax_trajectory(model, lag_order, spec, t, X(1:lag_order, :));
            else
                [t_sim, X_sim, sim_info] = simulate_trajectory(model, spec, [t(1), t(end)], X(1, :), t);
            end
            dyn_class = classify_dynamics(t_sim, X_sim, sim_info);
            sim_period = estimate_trajectory_period(t_sim, X_sim);
            if isfinite(sim_period.mean_period) && isfinite(real_period.mean_period) && real_period.mean_period > 0
                period_error = abs(sim_period.mean_period - real_period.mean_period) / real_period.mean_period;
            else
                period_error = NaN; % model has no clear period (fixed point, etc.) -- not comparable
            end
        else
            dyn_class = 'diverged';
            period_error = Inf; % diverged -- treat as maximally period-mismatched, not "not comparable"
        end

        results(end+1) = struct(... %#ok<AGROW>
            'variant', variant, 'flavor', run_c.flavor, 'lambda', lambda, ...
            'poly_order', poly_order_field, 'window_points', window_points, ...
            'test_function_order', test_function_order, 'lag_order', lag_order, ...
            'hill_n', hill_n_field, 'hill_K', hill_K_field, ...
            'library_spec', spec, 'num_active_terms', num_active, 'trajectory_rmse', metrics.rmse, ...
            'normalized_rmse', metrics.normalized_rmse, 'dynamics_class', dyn_class, ...
            'period_error', period_error, 'Xi', model, 'library_names', {library_names}, 'success', metrics.success);

    catch ME
        results(end+1) = struct(... %#ok<AGROW>
            'variant', variant, 'flavor', run_c.flavor, 'lambda', lambda, ...
            'poly_order', poly_order_field, 'window_points', window_points, ...
            'test_function_order', test_function_order, 'lag_order', lag_order, ...
            'hill_n', hill_n_field, 'hill_K', hill_K_field, ...
            'library_spec', spec, 'num_active_terms', NaN, 'trajectory_rmse', Inf, ...
            'normalized_rmse', Inf, 'dynamics_class', 'error', ...
            'period_error', Inf, 'Xi', [], 'library_names', {{}}, 'success', false);
        fprintf('  [WARN] combo %d/%d failed: %s\n', c, nCombos, ME.message);
    end

    if mod(c, progressStep) == 0
        fprintf('  ... %d/%d combinations complete\n', c, nCombos);
    end
end

fprintf('[GRID SEARCH] Done. %d/%d combinations succeeded.\n', sum([results.success]), nCombos);

end


function v = local_getfield_default(s, fname, default)
    if isfield(s, fname)
        v = s.(fname);
    else
        v = default;
    end
end


function [value, remaining] = local_extract_option(varargin_in, name, default)
    value = default;
    remaining = varargin_in;
    idx = find(strcmpi(remaining, name));
    if ~isempty(idx)
        value = remaining{idx(1) + 1};
        remaining([idx(1), idx(1) + 1]) = [];
    end
end