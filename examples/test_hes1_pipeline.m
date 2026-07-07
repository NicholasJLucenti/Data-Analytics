%% Test Script: Full Automated SINDy Pipeline on Hes1 Data
% Exercises the whole chain: io -> diagnostics -> selection ->
% preprocessing -> variant (dynamically routed) -> ensembling -> report,
% using the Hes1 gene regulatory network as the validation dataset.
clear; clc; close all;

%% 0. Configuration -- point this at wherever your Hes1 .mat files live
dataDir = fullfile(pwd, 'data');

mrnaValueFile = fullfile(dataDir, 'Hes1mRNAData.mat');
mrnaTimeFile  = fullfile(dataDir, 'Hes1mRNATime.mat');
protValueFile = fullfile(dataDir, 'Hes1Data.mat');
protTimeFile  = fullfile(dataDir, 'Hes1Time.mat');

lambda_standard = 0.35;
lambda_weak     = 0.35;  % weak-form system has a different natural scale than
                          % the pointwise system -- tune this independently,
                          % don't assume it should match lambda_standard
poly_order = 2;
ensemble_replicates = 200;

%% 1. Load raw data into the standard struct-array format
fprintf('Loading Hes1 datasets...\n');

mRNA_val = load(mrnaValueFile);
mRNA_t   = load(mrnaTimeFile);
prot_val = load(protValueFile);
prot_t   = load(protTimeFile);

raw(1).name = 'mRNA';
raw(1).t = mRNA_t.(char(fieldnames(mRNA_t)));
raw(1).y = mRNA_val.(char(fieldnames(mRNA_val)));
raw(1).source = mrnaValueFile;

raw(2).name = 'protein';
raw(2).t = prot_t.(char(fieldnames(prot_t)));
raw(2).y = prot_val.(char(fieldnames(prot_val)));
raw(2).source = protValueFile;

%% 2. Validate structural integrity before doing any real work
report = validate_input(raw);
if ~report.ok
    error('Input validation failed:\n%s', strjoin(report.messages, '\n'));
elseif ~isempty(report.messages)
    fprintf('Validation warnings:\n');
    for i = 1:numel(report.messages)
        fprintf('  %s\n', report.messages{i});
    end
end

%% 3. Diagnose each channel BEFORE merging (noise, sampling, sparsity, horizon)
fprintf('\nRunning diagnostics...\n');
diagnostics = profile_data(raw);
for i = 1:numel(diagnostics.per_variable)
    v = diagnostics.per_variable(i);
    fprintf('  -> %s: SNR = %.2f dB | dt_cv = %.4f | points/period = %.1f\n', ...
        v.name, v.noise.snr_db, v.sampling.dt_cv, v.sparsity.points_per_period);
end

%% 4. Decide preprocessing/variant route from those diagnostics
route = select_route(diagnostics);
fprintf('\n[ROUTE] Recommended variant: %s (%s)\n', route.recommended_variant, route.differentiation);
fprintf('        Reason: %s\n', route.reason);
for i = 1:numel(route.warnings)
    fprintf('        [WARNING] %s\n', route.warnings{i});
end

%% 5. Align, truncate, and densify onto one shared (t, X) grid
fprintf('\nAligning and densifying...\n');
target_density = 300;
[t, X, names] = align_and_truncate(raw, target_density);

%% 6. Run the recommended variant
% Standard SINDy needs smoothed data and finite-difference derivatives.
% Weak SINDy deliberately skips both -- it integrates the raw densified
% signal directly, which is the entire point of using it on noisy/sparse
% data. Do NOT smooth before calling run_weak_sindy.
fprintf('\nRunning SINDy solver (%s variant)...\n', route.recommended_variant);

if strcmp(route.recommended_variant, 'weak')
    [Xi, term_names, L, Y, weak_info] = run_weak_sindy(t, X, lambda_weak, poly_order);
    fprintf('  weak-form system: %d test functions, window = %d points, radius = %.4g\n', ...
        weak_info.num_test_functions, weak_info.window_points, weak_info.radius);
else
    X_smooth = smooth_data(t, X, 0.2);
    dXdt = compute_derivatives(t, X_smooth);
    [Xi, term_names, L, Y] = run_standard_sindy(X_smooth, dXdt, lambda_standard, poly_order);
end

%% 7. Ensemble: bootstrap the same (L, Y) system to get coefficient uncertainty
fprintf('\nRunning ensemble (%d bootstrap replicates)...\n', ensemble_replicates);
if strcmp(route.recommended_variant, 'weak')
    lambda_for_ensemble = lambda_weak;
else
    lambda_for_ensemble = lambda_standard;
end
ensemble = ensemble_sindy(L, Y, lambda_for_ensemble, ...
    'NumReplicates', ensemble_replicates, 'InclusionThreshold', 0.5);

%% 8. Report: single-shot fit vs. ensemble-aggregated fit, with inclusion probabilities
fprintf('\n=== SINGLE-SHOT DISCOVERED EQUATIONS (%s variant) ===\n', route.recommended_variant);
print_equations(Xi, term_names, names);

fprintf('\n=== ENSEMBLE-AGGREGATED EQUATIONS (median of stable terms, inclusion >= 50%%) ===\n');
print_equations(ensemble.Xi_final, term_names, names);

fprintf('\n--- Term inclusion probabilities across %d bootstrap replicates ---\n', ensemble.num_replicates);
for d = 1:size(X,2)
    fprintf('  %s:\n', names{d});
    for m = 1:numel(term_names)
        if ensemble.inclusion_prob(m, d) > 0
            fprintf('    %-10s  inclusion = %5.1f%%   mean = %+8.4f   std = %7.4f\n', ...
                term_names{m}, 100*ensemble.inclusion_prob(m, d), ...
                ensemble.Xi_mean(m, d), ensemble.Xi_std(m, d));
        end
    end
end

%% 9. Plot final clean trajectories
figure('Name', 'Hes1 Pipeline Output');
for d = 1:size(X,2)
    subplot(size(X,2), 1, d);
    plot(t, X(:,d), 'LineWidth', 2); ylabel(names{d}); grid on;
end
xlabel('Time');


%% --- local helper ---
function print_equations(Xi, term_names, state_names)
    for d = 1:size(Xi, 2)
        fprintf('d%s/dt = ', state_names{d});
        active_idx = find(Xi(:, d) ~= 0);
        if isempty(active_idx)
            fprintf('0\n');
        else
            for k = 1:length(active_idx)
                idx = active_idx(k);
                fprintf('%+.4f*%s ', Xi(idx, d), term_names{idx});
            end
            fprintf('\n');
        end
    end
end