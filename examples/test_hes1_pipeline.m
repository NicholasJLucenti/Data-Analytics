%% Test Script: Diagnostics on Local Mat-File Hes1 Data
clear; clc; close all;

% Add repository folders to MATLAB path
addpath('C:\Users\nickj\MATLAB Drive\Data Analytics\Data-Analytics-main\diagnostics');


fprintf('Loading Hes1 .mat datasets...\n');

%% 1. Load the Data
% MATLAB loads the variables directly into the workspace. 
% (Note: If your internal variable names differ from the file names, 
% change the right side of the assignments below accordingly.)
mRNA_struct   = load("C:\Users\nickj\MATLAB Drive\Compiled Works\Hes1 Data\Hes1mRNAData.mat");
mRNA_t_struct = load("C:\Users\nickj\MATLAB Drive\Compiled Works\Hes1 Data\Hes1mRNATime.mat");
prot_struct   = load("C:\Users\nickj\MATLAB Drive\Compiled Works\Hes1 Data\Hes1Data.mat");
prot_t_struct = load("C:\Users\nickj\MATLAB Drive\Compiled Works\Hes1 Data\Hes1Time.mat");


mrna_raw   = mRNA_struct.(char(fields(mRNA_struct)));
t_mrna_raw = mRNA_t_struct.(char(fields(mRNA_t_struct)));
prot_raw   = prot_struct.(char(fields(prot_struct)));
t_prot_raw = prot_t_struct.(char(fields(prot_t_struct)));

%% Phase 1 & 2: Process, Align, Truncate, and Densify Data Automatically
fprintf('\nExecuting Preprocessing Layer...\n');
target_density = 300; % Bumping up the data volume to 300 points
[t, X] = align_and_truncate(t_mrna_raw, mrna_raw, t_prot_raw, prot_raw, target_density);

%% Phase 3: Run System Quality Profiler (Will now see 300 samples!)
report = profile_data(t, X);

%% Phase 4: Compute Governing Equations via SINDy Solver
fprintf('\nExecuting SINDy Solver Layer...\n');
lambda = 0.35;      % Raised lambda slightly to ensure sparsity on smooth trends
poly_order = 2;    

[Xi, term_names] = run_standard_sindy(t, X, lambda, poly_order);

% Print Discovered Equations
fprintf('\n=== DISCOVERED GOVERNING EQUATIONS ===\n');
for d = 1:size(X,2)
    fprintf('dx%d/dt = ', d);
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

%% Plot Final Clean Visualizations
figure('Name', 'Hes1 Perfect Alignment');
subplot(2,1,1); plot(t, X(:,1), 'b-', 'LineWidth', 2); ylabel('mRNA'); grid on;
subplot(2,1,2); plot(t, X(:,2), 'r-', 'LineWidth', 2); ylabel('Protein'); xlabel('Time'); grid on;
