function best = select_best_model(results, varargin)
%SELECT_BEST_MODEL Extract representative models from a grid-search
%results table along several distinct, sometimes-competing criteria.
%
%   best = SELECT_BEST_MODEL(results, ...)
%
%   Input:
%     results - struct array from benchmarking/run_grid_search.m
%
%   Name-value options:
%     'ErrorTolerance'     - for the "sparsest" pick: how much worse
%                             (relative) trajectory RMSE than the single
%                             best RMSE in the whole sweep is still
%                             acceptable (default 0.25, i.e. within 25%
%                             of the best RMSE). Without this, "sparsest"
%                             would trivially pick a near-empty model that
%                             fits nothing.
%     'RequireOscillatory' - if true (default), the "topology preserving"
%                             pick is restricted to
%                             dynamics_class == 'oscillatory'. Set false
%                             if your dataset's real dynamics are a fixed
%                             point rather than a limit cycle -- in that
%                             case you likely want 'fixed_point' instead;
%                             this option only controls the built-in
%                             oscillatory case, see note below for others.
%
%   Output: best struct with fields:
%     .lowest_error        - result with minimum trajectory_rmse among
%                             successful runs
%     .sparsest             - fewest num_active_terms among results
%                             within ErrorTolerance of the best RMSE
%     .topology_preserving  - lowest RMSE among results classified
%                             'oscillatory' ([] if none found)
%     .all_valid            - the filtered results actually used
%                             (success == true, finite RMSE), for any
%                             further custom analysis (e.g. picking
%                             'fixed_point' models by hand via
%                             best.all_valid(strcmp({best.all_valid.dynamics_class},'fixed_point')))

p = inputParser;
addParameter(p, 'ErrorTolerance', 0.25, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'RequireOscillatory', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

successMask = [results.success] & isfinite([results.trajectory_rmse]);
valid = results(successMask);

if isempty(valid)
    error('select_best_model:noValidModels', ...
        'No grid-search combination produced a finite, successful trajectory fit.');
end

rmseVals = [valid.trajectory_rmse];
[bestRmse, bestIdx] = min(rmseVals);
best.lowest_error = valid(bestIdx);

sparseCandidates = valid(rmseVals <= bestRmse * (1 + opts.ErrorTolerance));
[~, sparseLocalIdx] = min([sparseCandidates.num_active_terms]);
best.sparsest = sparseCandidates(sparseLocalIdx);

if opts.RequireOscillatory
    oscMask = strcmp({valid.dynamics_class}, 'oscillatory');
    oscCandidates = valid(oscMask);
    if isempty(oscCandidates)
        best.topology_preserving = [];
        warning('select_best_model:noOscillatoryModel', ...
            'No swept model sustained oscillatory dynamics under forward simulation.');
    else
        oscRmse = [oscCandidates.trajectory_rmse];
        [~, oscLocalIdx] = min(oscRmse);
        best.topology_preserving = oscCandidates(oscLocalIdx);
    end
else
    best.topology_preserving = [];
end

best.all_valid = valid;

end