%% Test Script: Full Automated SINDy Pipeline on Hes1 Data
% Exercises the whole chain: io -> diagnostics -> selection ->
% preprocessing -> variant -> report, using the Hes1 gene regulatory
% network as the validation dataset.
clear; clc; close all;

%% 0. Configuration -- point this at wherever your Hes1 .mat files live
dataDir = fullfile(pwd, 'data');

mrnaValueFile = fullfile(dataDir, 'Hes1mRNAData.mat');
mrnaTimeFile  = fullfile(dataDir, 'Hes1mRNATime.mat');
protValueFile = fullfile(dataDir, 'Hes1Data.mat');
protTimeFile  = fullfile(dataDir, 'Hes1Time.mat');

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
if ~route.implemented
    fprintf(['        NOTE: "%s" is not implemented in variants/ yet. ' ...
        'Falling back to run_standard_sindy.m for this run.\n'], route.recommended_variant);
end

%% 5. Align, truncate, and densify onto one shared (t, X) grid
fprintf('\nAligning and densifying...\n');
target_density = 300;
[t, X, names] = align_and_truncate(raw, target_density);

%% 6. Shared preprocessing: smooth once, differentiate once
X_smooth = smooth_data(t, X, 0.2);
dXdt = compute_derivatives(t, X_smooth);

%% 7. Run SINDy (standard variant, until weak/implicit exist)
fprintf('\nRunning SINDy solver...\n');
lambda = 0.35;
poly_order = 2;
[Xi, term_names] = run_standard_sindy(X_smooth, dXdt, lambda, poly_order);

%% 8. Report discovered equations
fprintf('\n=== DISCOVERED GOVERNING EQUATIONS ===\n');
for d = 1:size(X,2)
    fprintf('d%s/dt = ', names{d});
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

%% 9. Plot final clean trajectories
figure('Name', 'Hes1 Pipeline Output');
for d = 1:size(X,2)
    subplot(size(X,2), 1, d);
    plot(t, X(:,d), 'LineWidth', 2); ylabel(names{d}); grid on;
end
xlabel('Time');