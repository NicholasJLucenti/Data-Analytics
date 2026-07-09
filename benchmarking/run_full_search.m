function [results, t, X, names] = run_full_search(raw, varargin)
%RUN_FULL_SEARCH Sweep every implemented SINDy variant (standard, weak,
%implicit), each across its own library-flavor grid (plain polynomial and
%Hill/saturation flavors), and merge into one combined results table.
%
%   [results, t, X, names] = RUN_FULL_SEARCH(raw, ...)
%
%   Name-value options:
%     'Resolution' - 'fast' | 'balanced' (default) | 'thorough'. Sets ALL
%                     grid defaults at once via get_resolution_preset.m;
%                     any grid also passed explicitly overrides just that
%                     field. E.g. run_full_search(raw, 'Resolution', 'fast')
%                     for a quick pass, or 'thorough' once you trust the
%                     pipeline and want the best answer. See
%                     get_resolution_preset.m for approximate combination
%                     counts per level.
%     'LambdaGrid', 'PolyOrderGrid', 'WindowPointsGrid',
%     'TestFunctionOrderGrid', 'NumDensePoints', 'SmoothingFactor',
%     'LibraryFlavorGrid', 'HillCoeffGrid', 'HillKCandidates'
%         - forwarded to run_grid_search.m for each variant; override
%           individual fields of whichever Resolution preset is active
%     'Variants' - cell array of variant names to sweep
%                  (default {'standard', 'weak', 'implicit'})
%
%   Output:
%     results - merged struct array from every swept variant, each entry
%               tagged with .variant (see run_grid_search.m for full
%               field list)
%     t, X, names - the shared densified trajectory (identical across
%                   variants since they share raw/NumDensePoints)

[resolution, varargin] = local_extract_option(varargin, 'Resolution', 'balanced');
preset = get_resolution_preset(resolution);

p = inputParser;
addParameter(p, 'LambdaGrid', preset.LambdaGrid, @isnumeric);
addParameter(p, 'PolyOrderGrid', preset.PolyOrderGrid, @isnumeric);
addParameter(p, 'WindowPointsGrid', preset.WindowPointsGrid, @isnumeric);
addParameter(p, 'TestFunctionOrderGrid', preset.TestFunctionOrderGrid, @isnumeric);
addParameter(p, 'NumDensePoints', preset.NumDensePoints, @isnumeric);
addParameter(p, 'SmoothingFactor', 0.2, @isnumeric);
addParameter(p, 'LibraryFlavorGrid', preset.LibraryFlavorGrid, @iscell);
addParameter(p, 'HillCoeffGrid', preset.HillCoeffGrid, @isnumeric);
addParameter(p, 'HillKCandidates', preset.HillKCandidates, @isnumeric);
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

fprintf('[FULL SEARCH] Resolution: %s\n', resolution);

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


function [value, remaining] = local_extract_option(varargin_in, name, default)
    value = default;
    remaining = varargin_in;
    idx = find(strcmpi(remaining, name));
    if ~isempty(idx)
        value = remaining{idx(1) + 1};
        remaining([idx(1), idx(1) + 1]) = [];
    end
end