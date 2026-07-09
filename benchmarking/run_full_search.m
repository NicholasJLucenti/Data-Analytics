function [results, t, X, names] = run_full_search(raw, varargin)
%RUN_FULL_SEARCH Sweep every implemented SINDy variant (standard, weak,
%implicit) and merge into one combined results table, so
%select_best_model.m picks the best model across variants, not just
%within one.
%
%   [results, t, X, names] = RUN_FULL_SEARCH(raw, ...)
%
%   Input:
%     raw - raw struct array (io/load_raw_data.m format)
%
%   Name-value options (forwarded to run_grid_search.m for each variant):
%     'LambdaGrid', 'PolyOrderGrid', 'WindowPointsGrid',
%     'TestFunctionOrderGrid', 'NumDensePoints', 'SmoothingFactor'
%     'Variants' - cell array of variant names to sweep
%                  (default {'standard', 'weak', 'implicit'})
%
%   Output:
%     results - merged struct array from every swept variant, each entry
%               tagged with .variant (see run_grid_search.m for full
%               field list)
%     t, X, names - the shared densified trajectory (identical across
%                   variants since they share raw/NumDensePoints)
%
%   Runtime note: sweeping all three variants multiplies total combinations.
%   With default grids that's roughly (|Lambda|*|Poly|) for standard,
%   (|Lambda|*|Poly|) for implicit, and (|Lambda|*|Poly|*|Window|*|TFO|)
%   for weak -- weak dominates the total runtime. Shrink LambdaGrid or the
%   weak-specific grids first if this is too slow.

p = inputParser;
addParameter(p, 'LambdaGrid', logspace(-2, 0, 8), @isnumeric);
addParameter(p, 'PolyOrderGrid', [1 2 3], @isnumeric);
addParameter(p, 'WindowPointsGrid', [11 21 31], @isnumeric);
addParameter(p, 'TestFunctionOrderGrid', [2 4 6], @isnumeric);
addParameter(p, 'NumDensePoints', 300, @isnumeric);
addParameter(p, 'SmoothingFactor', 0.2, @isnumeric);
addParameter(p, 'Variants', {'standard', 'weak', 'implicit'}, @iscell);
parse(p, varargin{:});
opts = p.Results;

results = struct('variant', {}, 'lambda', {}, 'poly_order', {}, 'window_points', {}, 'test_function_order', {}, ...
    'num_active_terms', {}, 'trajectory_rmse', {}, 'normalized_rmse', {}, 'dynamics_class', {}, ...
    'Xi', {}, 'library_names', {}, 'success', {});

t = [];
X = [];
names = {};

for i = 1:numel(opts.Variants)
    variant = opts.Variants{i};
    fprintf('\n[FULL SEARCH] === Variant: %s ===\n', variant);
    [r, t, X, names] = run_grid_search(raw, variant, ...
        'LambdaGrid', opts.LambdaGrid, 'PolyOrderGrid', opts.PolyOrderGrid, ...
        'WindowPointsGrid', opts.WindowPointsGrid, 'TestFunctionOrderGrid', opts.TestFunctionOrderGrid, ...
        'NumDensePoints', opts.NumDensePoints, 'SmoothingFactor', opts.SmoothingFactor);
    results = [results, r]; %#ok<AGROW>
end

fprintf('\n[FULL SEARCH] Done. %d total combinations across %d variant(s) (%d succeeded).\n', ...
    numel(results), numel(opts.Variants), sum([results.success]));

end