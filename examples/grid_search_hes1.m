%% Grid Search: Sweep SINDy hyperparameters on Hes1 data and extract
%% competing best-fit models (sparsest, lowest error, topology-preserving).
clear; clc; close all;

%% 0. Configuration
dataDir = fullfile(pwd, 'data');
mrnaValueFile = fullfile(dataDir, 'Hes1mRNAData.mat');
mrnaTimeFile  = fullfile(dataDir, 'Hes1mRNATime.mat');
protValueFile = fullfile(dataDir, 'Hes1Data.mat');
protTimeFile  = fullfile(dataDir, 'Hes1Time.mat');

%% 1. Load raw data
mRNA_val = load(mrnaValueFile);
mRNA_t   = load(mrnaTimeFile);
prot_val = load(protValueFile);
prot_t   = load(protTimeFile);

raw(1).name = 'mRNA';
raw(1).t = mRNA_t.(char(fieldnames(mRNA_t)));
raw(1).y = mRNA_val.(char(fieldnames(mRNA_val)));

raw(2).name = 'protein';
raw(2).t = prot_t.(char(fieldnames(prot_t)));
raw(2).y = prot_val.(char(fieldnames(prot_val)));

%% 2. Diagnose + route (as before) to decide which variant to sweep
diagnostics = profile_data(raw);
route = select_route(diagnostics);
fprintf('[ROUTE] Recommended variant: %s\n', route.recommended_variant);

%% 3. Run the full grid search on the recommended variant
% This is the "most thorough" sweep: lambda x poly_order x, for weak-SINDy,
% window_points x test_function_order. This is the slow, expensive step --
% expect this to take a while since every combination refits from scratch
% and forward-simulates for scoring.
[results, t, X, names] = run_grid_search(raw, route.recommended_variant, ...
    'LambdaGrid', logspace(-2, 0, 8), ...
    'PolyOrderGrid', [1 2 3], ...
    'WindowPointsGrid', [11 21 31], ...
    'TestFunctionOrderGrid', [2 4 6]);

%% 4. Extract competing best models
best = select_best_model(results, 'ErrorTolerance', 0.25, 'RequireOscillatory', true);

fprintf('\n=== LOWEST TRAJECTORY ERROR MODEL ===\n');
report_model(best.lowest_error, names);

fprintf('\n=== SPARSEST MODEL (within 25%% of best RMSE) ===\n');
report_model(best.sparsest, names);

if ~isempty(best.topology_preserving)
    fprintf('\n=== TOPOLOGY-PRESERVING MODEL (sustains oscillation) ===\n');
    report_model(best.topology_preserving, names);
else
    fprintf('\n=== TOPOLOGY-PRESERVING MODEL ===\n');
    fprintf('No swept combination sustained oscillatory dynamics under forward simulation.\n');
    fprintf('Consider widening the grid, or check whether the real data itself is actually\n');
    fprintf('oscillatory (see diagnostics.per_variable(...).sparsity) before trusting any\n');
    fprintf('single-shot fit''s qualitative behavior.\n');
end

%% 5. Quick sparsity-vs-error trade-off view across the whole valid sweep
figure('Name', 'Grid Search: Sparsity vs. Trajectory Error');
scatter([best.all_valid.num_active_terms], [best.all_valid.trajectory_rmse], 20, 'filled');
xlabel('Number of active terms (sparsity)');
ylabel('Trajectory RMSE');
title('Grid search results: complexity vs. fit trade-off');
set(gca, 'YScale', 'log');
grid on;

%% 6. Simulate the extracted models forward and overlay against real data
plot_best_models(t, X, names, best);


%% --- local helper ---
function report_model(r, names)
    fprintf('  lambda = %.4g | poly_order = %d', r.lambda, r.poly_order);
    if ~isnan(r.window_points)
        fprintf(' | window_points = %d | test_function_order = %d', r.window_points, r.test_function_order);
    end
    fprintf('\n  active terms = %d | trajectory RMSE = %.4f | dynamics = %s\n', ...
        r.num_active_terms, r.trajectory_rmse, r.dynamics_class);
    for d = 1:size(r.Xi, 2)
        fprintf('  d%s/dt = ', names{d});
        active_idx = find(r.Xi(:, d) ~= 0);
        if isempty(active_idx)
            fprintf('0\n');
        else
            for k = 1:length(active_idx)
                idx = active_idx(k);
                fprintf('%+.4f*%s ', r.Xi(idx, d), r.library_names{idx});
            end
            fprintf('\n');
        end
    end
end