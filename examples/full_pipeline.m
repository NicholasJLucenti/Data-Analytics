%% Full Automated Pipeline: multi-variant grid search + cross-variant selection
% Single entry point: load raw data, diagnose it, sweep standard/weak/
% implicit SINDy together, then extract and plot the best models across
% ALL variants at once (not just one fixed choice).
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

%% 2. Validate
report = validate_input(raw);
if ~report.ok
    error('Input validation failed:\n%s', strjoin(report.messages, '\n'));
elseif ~isempty(report.messages)
    fprintf('Validation warnings:\n');
    for i = 1:numel(report.messages)
        fprintf('  %s\n', report.messages{i});
    end
end

%% 3. Diagnose (informational -- full search below tries every variant
%% regardless, so this is no longer a hard gate, just context)
diagnostics = profile_data(raw);
for i = 1:numel(diagnostics.per_variable)
    v = diagnostics.per_variable(i);
    fprintf('  -> %s: SNR = %.2f dB | dt_cv = %.4f | points/period = %.1f\n', ...
        v.name, v.noise.snr_db, v.sampling.dt_cv, v.sparsity.points_per_period);
end
route = select_route(diagnostics);
fprintf('[INFO] Single-variant diagnostic recommendation: %s (informational -- full search covers all variants)\n', ...
    route.recommended_variant);

%% 4. Sweep ALL implemented variants together
[results, t, X, names] = run_full_search(raw, ...
    'LambdaGrid', logspace(-2, 0, 6), ...
    'PolyOrderGrid', [1 2 3], ...
    'WindowPointsGrid', [11 21 31], ...
    'TestFunctionOrderGrid', [2 4 6]);

%% 5. Extract best models across ALL variants
best = select_best_model(results, 'ErrorTolerance', 0.25, 'RequireOscillatory', true);

fprintf('\n=== LOWEST TRAJECTORY ERROR MODEL (variant: %s) ===\n', best.lowest_error.variant);
print_equations(best.lowest_error.Xi, best.lowest_error.library_names, names);
fprintf('  active terms = %d | rmse = %.4f | dynamics = %s\n', ...
    best.lowest_error.num_active_terms, best.lowest_error.trajectory_rmse, best.lowest_error.dynamics_class);

fprintf('\n=== SPARSEST MODEL (variant: %s) ===\n', best.sparsest.variant);
print_equations(best.sparsest.Xi, best.sparsest.library_names, names);
fprintf('  active terms = %d | rmse = %.4f | dynamics = %s\n', ...
    best.sparsest.num_active_terms, best.sparsest.trajectory_rmse, best.sparsest.dynamics_class);

if ~isempty(best.topology_preserving)
    fprintf('\n=== TOPOLOGY-PRESERVING MODEL (variant: %s) ===\n', best.topology_preserving.variant);
    print_equations(best.topology_preserving.Xi, best.topology_preserving.library_names, names);
    fprintf('  active terms = %d | rmse = %.4f | dynamics = %s\n', ...
        best.topology_preserving.num_active_terms, best.topology_preserving.trajectory_rmse, best.topology_preserving.dynamics_class);
else
    fprintf('\nNo swept model (any variant) sustained oscillatory dynamics under forward simulation.\n');
end

%% 6. Sparsity-vs-error trade-off across the whole valid sweep, colored by variant
figure('Name', 'Full Search: Sparsity vs. Trajectory Error by Variant');
hold on;
variantColors = containers.Map({'standard', 'weak', 'implicit'}, {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]});
variantList = unique({best.all_valid.variant});
for i = 1:numel(variantList)
    vMask = strcmp({best.all_valid.variant}, variantList{i});
    subset = best.all_valid(vMask);
    scatter([subset.num_active_terms], [subset.trajectory_rmse], 20, ...
        variantColors(variantList{i}), 'filled', 'DisplayName', variantList{i});
end
xlabel('Number of active terms (sparsity)');
ylabel('Trajectory RMSE');
title('Full search results: complexity vs. fit, by variant');
set(gca, 'YScale', 'log');
legend('Location', 'best');
grid on;

%% 7. Simulate the extracted models forward and overlay against real data
plot_best_models(t, X, names, best);