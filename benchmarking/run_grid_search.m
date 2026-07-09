function [results, t, X, names] = run_grid_search(raw, variant, varargin)
%RUN_GRID_SEARCH Sweep SINDy hyperparameters and score every resulting
%model, so a downstream selection step (select_best_model.m) can pick
%among competing trade-offs -- sparsest, lowest trajectory error,
%topology-preserving -- instead of committing to one fixed-parameter fit.
%
%   [results, t, X, names] = RUN_GRID_SEARCH(raw, variant, ...)
%
%   Inputs:
%     raw     - raw struct array (io/load_raw_data.m format)
%     variant - 'standard' or 'weak'
%
%   Name-value options:
%     'LambdaGrid'            - sparsity thresholds to sweep
%                                (default logspace(-2, 0, 8))
%     'PolyOrderGrid'         - candidate library orders to sweep
%                                (default [1 2 3])
%     'WindowPointsGrid'      - (weak only) test-function window sizes
%                                (default [11 21 31])
%     'TestFunctionOrderGrid' - (weak only) test-function smoothness
%                                orders (default [2 4 6])
%     'NumDensePoints'        - densification target for
%                                align_and_truncate (default 300)
%     'SmoothingFactor'       - (standard only) Lowess span fraction
%                                (default 0.2)
%
%   Output:
%     results - struct array, one entry per parameter combination:
%                 .lambda, .poly_order, .window_points, .test_function_order
%                   (window_points/test_function_order are NaN for the
%                   standard variant, which doesn't use them)
%                 .num_active_terms  - nnz(Xi), a sparsity count
%                 .trajectory_rmse   - from compute_trajectory_error.m
%                 .normalized_rmse   - scale-free version of the above
%                 .dynamics_class    - from classify_dynamics.m
%                 .Xi, .library_names - the fitted model itself
%                 .success           - false if this combination errored
%                                       or failed to simulate
%     t, X, names - the shared densified trajectory/channel names used to
%                   score every candidate (from align_and_truncate.m),
%                   returned so callers don't have to recompute them

p = inputParser;
addParameter(p, 'LambdaGrid', logspace(-2, 0, 8), @isnumeric);
addParameter(p, 'PolyOrderGrid', [1 2 3], @isnumeric);
addParameter(p, 'WindowPointsGrid', [11 21 31], @isnumeric);
addParameter(p, 'TestFunctionOrderGrid', [2 4 6], @isnumeric);
addParameter(p, 'NumDensePoints', 300, @isnumeric);
addParameter(p, 'SmoothingFactor', 0.2, @isnumeric);
parse(p, varargin{:});
opts = p.Results;

if ~ismember(variant, {'standard', 'weak'})
    error('run_grid_search:badVariant', 'variant must be ''standard'' or ''weak''.');
end

%% Shared preprocessing -- done once, independent of the swept parameters
[t, X, names] = align_and_truncate(raw, opts.NumDensePoints);

if strcmp(variant, 'standard')
    X_smooth = smooth_data(t, X, opts.SmoothingFactor);
    dXdt = compute_derivatives(t, X_smooth);
end

%% Build the parameter grid
if strcmp(variant, 'standard')
    [lambdaGrid, polyGrid] = ndgrid(opts.LambdaGrid, opts.PolyOrderGrid);
    combos = [lambdaGrid(:), polyGrid(:), nan(numel(lambdaGrid), 2)];
else
    [lambdaGrid, polyGrid, winGrid, tfoGrid] = ndgrid(opts.LambdaGrid, opts.PolyOrderGrid, ...
        opts.WindowPointsGrid, opts.TestFunctionOrderGrid);
    combos = [lambdaGrid(:), polyGrid(:), winGrid(:), tfoGrid(:)];
end
nCombos = size(combos, 1);

fprintf('[GRID SEARCH] Sweeping %d parameter combinations (%s variant)...\n', nCombos, variant);

results = struct('lambda', {}, 'poly_order', {}, 'window_points', {}, 'test_function_order', {}, ...
    'num_active_terms', {}, 'trajectory_rmse', {}, 'normalized_rmse', {}, 'dynamics_class', {}, ...
    'Xi', {}, 'library_names', {}, 'success', {});

progressStep = max(1, round(nCombos / 10));

for c = 1:nCombos
    lambda = combos(c, 1);
    poly_order = combos(c, 2);
    window_points = combos(c, 3);
    test_function_order = combos(c, 4);

    try
        if strcmp(variant, 'standard')
            [Xi, library_names] = run_standard_sindy(X_smooth, dXdt, lambda, poly_order);
        else
            [Xi, library_names] = run_weak_sindy(t, X, lambda, poly_order, ...
                'WindowPoints', window_points, 'TestFunctionOrder', test_function_order);
        end

        metrics = compute_trajectory_error(t, X, Xi, poly_order);

        if metrics.success
            [t_sim, X_sim, sim_info] = simulate_trajectory(Xi, poly_order, [t(1), t(end)], X(1,:), t);
            dyn_class = classify_dynamics(t_sim, X_sim, sim_info);
        else
            dyn_class = 'diverged';
        end

        results(end+1) = struct(... %#ok<AGROW>
            'lambda', lambda, 'poly_order', poly_order, ...
            'window_points', window_points, 'test_function_order', test_function_order, ...
            'num_active_terms', nnz(Xi), 'trajectory_rmse', metrics.rmse, ...
            'normalized_rmse', metrics.normalized_rmse, 'dynamics_class', dyn_class, ...
            'Xi', Xi, 'library_names', {library_names}, 'success', metrics.success);

    catch ME
        results(end+1) = struct(... %#ok<AGROW>
            'lambda', lambda, 'poly_order', poly_order, ...
            'window_points', window_points, 'test_function_order', test_function_order, ...
            'num_active_terms', NaN, 'trajectory_rmse', Inf, ...
            'normalized_rmse', Inf, 'dynamics_class', 'error', ...
            'Xi', [], 'library_names', {{}}, 'success', false);
        fprintf('  [WARN] combo %d/%d failed: %s\n', c, nCombos, ME.message);
    end

    if mod(c, progressStep) == 0
        fprintf('  ... %d/%d combinations complete\n', c, nCombos);
    end
end

fprintf('[GRID SEARCH] Done. %d/%d combinations succeeded.\n', sum([results.success]), nCombos);

end