function fig = plot_best_models(t, X, names, best, varargin)
%PLOT_BEST_MODELS Forward-simulate the models extracted by
%select_best_model.m and overlay them against the real densified
%trajectory, one subplot per state channel.
%
%   fig = PLOT_BEST_MODELS(t, X, names, best, ...)
%
%   Inputs:
%     t, X, names - the same densified trajectory the grid search scored
%                   against (from preprocessing/align_and_truncate.m, or
%                   the 2nd/3rd/4th outputs of benchmarking/run_grid_search.m)
%     best        - struct from benchmarking/select_best_model.m, with
%                   fields .lowest_error, .sparsest, .topology_preserving
%                   (any of which may be empty and will be skipped)
%
%   Name-value options:
%     'MaxWallSeconds' - passed through to simulate_trajectory.m (default
%                        5 -- more generous than the grid-search default
%                        of 3, since this is a one-off plot, not one call
%                        among hundreds)
%
%   If two of the extracted models happen to be the identical parameter
%   combination (common -- e.g. "lowest error" and "topology preserving"
%   are often literally the same winning model), they are merged into a
%   single plotted line with a combined label instead of drawing two
%   identical overlapping curves.
%
%   A model that diverges during simulation is still drawn (dotted,
%   truncated at the point it was cut off) rather than silently dropped,
%   so a bad candidate is visually obvious as "runs off the trajectory"
%   rather than just missing.

p = inputParser;
addParameter(p, 'MaxWallSeconds', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});
opts = p.Results;

D = size(X, 2);

colorMap.lowest_error = [0.85, 0.33, 0.10];
colorMap.sparsest = [0.00, 0.45, 0.74];
colorMap.topology_preserving = [0.47, 0.67, 0.19];

candidateDefs = {
    'lowest_error',        'Lowest error';
    'sparsest',             'Sparsest';
    'topology_preserving',  'Topology-preserving'
};

keyList = {};
candidates = struct('label', {}, 'result', {}, 'color', {});

for i = 1:size(candidateDefs, 1)
    fn = candidateDefs{i, 1};
    label = candidateDefs{i, 2};
    if ~isfield(best, fn) || isempty(best.(fn))
        continue
    end
    r = best.(fn);
    key = sprintf('%.6g_%d_%.6g_%.6g', r.lambda, r.poly_order, r.window_points, r.test_function_order);
    matchIdx = find(strcmp(keyList, key), 1);
    if ~isempty(matchIdx)
        candidates(matchIdx).label = [candidates(matchIdx).label ' & ' label];
    else
        keyList{end+1} = key; %#ok<AGROW>
        candidates(end+1).label = label; %#ok<AGROW>
        candidates(end).result = r;
        candidates(end).color = colorMap.(fn);
    end
end

if isempty(candidates)
    error('plot_best_models:noModels', 'best struct contains no models to plot.');
end

fig = figure('Name', 'Model Fit Comparison');
for d = 1:D
    subplot(D, 1, d);
    hold on;
    plot(t, X(:, d), 'k-', 'LineWidth', 2.5, 'DisplayName', 'Real (densified) data');

    for c = 1:numel(candidates)
        r = candidates(c).result;
        [t_sim, X_sim, sim_info] = simulate_trajectory(r.Xi, r.library_spec, [t(1), t(end)], X(1, :), t, ...
            'MaxWallSeconds', opts.MaxWallSeconds);

        if sim_info.success && size(X_sim, 1) == numel(t)
            plot(t_sim, X_sim(:, d), '--', 'Color', candidates(c).color, 'LineWidth', 1.75, ...
                'DisplayName', sprintf('%s (RMSE=%.3f, %d terms)', candidates(c).label, r.trajectory_rmse, r.num_active_terms));
        elseif ~isempty(X_sim)
            plot(t_sim, X_sim(:, d), ':', 'Color', candidates(c).color, 'LineWidth', 1.25, ...
                'DisplayName', sprintf('%s (diverged: %s)', candidates(c).label, sim_info.message));
        end
    end

    ylabel(names{d});
    grid on;
    if d == 1
        legend('Location', 'best', 'Interpreter', 'none');
    end
end
xlabel('Time');

end