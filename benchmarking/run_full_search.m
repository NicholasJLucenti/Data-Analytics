function [results, t, X, names] = run_full_search(raw, varargin)
%RUN_FULL_SEARCH Sweep every implemented SINDy variant (standard, weak,
%implicit), each across its own library-flavor grid (plain polynomial and
%Hill/saturation flavors), and merge into one combined results table.
%
%   [results, t, X, names] = RUN_FULL_SEARCH(raw, ...)
%
%   Name-value options (forwarded to run_grid_search.m for each variant):
%     'LambdaGrid', 'PolyOrderGrid', 'WindowPointsGrid',
%     'TestFunctionOrderGrid', 'NumDensePoints', 'SmoothingFactor',
%     'LibraryFlavorGrid', 'HillCoeffGrid', 'HillKCandidates'
%     'Variants' - cell array of variant names to sweep
%                  (default {'standard', 'weak', 'implicit'})
%
%   See run_grid_search.m for defaults and the runtime note on Hill
%   flavor combinatorics -- the same conservative defaults apply here,
%   multiplied across however many variants are swept.

p = inputParser;
addParameter(p, 'LambdaGrid', logspace(-2, 0, 5), @isnumeric);
addParameter(p, 'PolyOrderGrid', [1 2 3], @isnumeric);
addParameter(p, 'WindowPointsGrid', [15 25], @isnumeric);
addParameter(p, 'TestFunctionOrderGrid', [2 4], @isnumeric);
addParameter(p, 'NumDensePoints', 300, @isnumeric);
addParameter(p, 'SmoothingFactor', 0.2, @isnumeric);
addParameter(p, 'LibraryFlavorGrid', {'poly_only', 'hill_mixed'}, @iscell);
addParameter(p, 'HillCoeffGrid', [1 2 4], @isnumeric);
addParameter(p, 'HillKCandidates', 3, @isnumeric);
addParameter(p, 'Variants', {'standard', 'weak', 'implicit'}, @iscell);
parse(p, varargin{:});
opts = p.Results;

results = struct('variant', {}, 'flavor', {}, 'lambda', {}, 'poly_order', {}, 'window_points', {}, ...
    'test_function_order', {}, 'hill_n', {}, 'hill_K', {}, 'library_spec', {}, ...
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
        'NumDensePoints', opts.NumDensePoints, 'SmoothingFactor', opts.SmoothingFactor, ...
        'LibraryFlavorGrid', opts.LibraryFlavorGrid, 'HillCoeffGrid', opts.HillCoeffGrid, ...
        'HillKCandidates', opts.HillKCandidates);
    results = [results, r]; %#ok<AGROW>
end

fprintf('\n[FULL SEARCH] Done. %d total combinations across %d variant(s) (%d succeeded).\n', ...
    numel(results), numel(opts.Variants), sum([results.success]));

end