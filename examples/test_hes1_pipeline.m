%% Test Script: Diagnostics on Local Mat-File Hes1 Data
clear; clc; close all;

% Add repository folders to MATLAB path
addpath('../diagnostics');
addpath('../preprocessing');

fprintf('Loading Hes1 .mat datasets...\n');

%% 1. Load the Data
% MATLAB loads the variables directly into the workspace. 
% (Note: If your internal variable names differ from the file names, 
% change the right side of the assignments below accordingly.)
mRNA_struct   = load('../data/Hes1mRNAData.mat');
mRNA_t_struct = load('../data/Hes1mRNATime.mat');
prot_struct   = load('../data/Hes1Data.mat');
prot_t_struct = load('../data/Hes1Time.mat');

% Dynamically extract the first variable found in each mat file
mrna_fields = fields(mRNA_struct);
mrna_raw    = mRNA_struct.(mrna_fields{1});

mrna_t_fields = fields(mRNA_t_struct);
t_mrna_raw    = mRNA_t_struct.(mrna_t_fields{1});

prot_fields = fields(prot_struct);
prot_raw    = prot_struct.(prot_fields{1});

prot_t_fields = fields(prot_t_struct);
t_prot_raw    = prot_t_struct.(prot_t_fields{1});

%% 2. Handle Replicates (Multiple measurements at the same time)
% We find unique time points and compute the mean concentration for each.
fprintf('Conditioning data: Averaging concurrent replicate measurements...\n');

[t_mrna_uni, ~, idx_mrna] = unique(t_mrna_raw);
mrna_clean = accumarray(idx_mrna, mrna_raw, [], @mean);

[t_prot_uni, ~, idx_prot] = unique(t_prot_raw);
prot_clean = accumarray(idx_prot, prot_raw, [], @mean);

%% 3. Align the Arrays
% Ensure both variables are mapped to the exact same time grid.
% If their unique time vectors don't match perfectly, we interpolate to the mRNA time grid.
if ~isequal(t_mrna_uni, t_prot_uni)
    fprintf('[INFO] Time grids mismatch slightly. Aligning protein data to mRNA time grid...\n');
    prot_aligned = interp1(t_prot_uni, prot_clean, t_mrna_uni, 'spline');
    t = t_mrna_uni;
    X = [mrna_clean, prot_aligned];
else
    t = t_mrna_uni;
    X = [mrna_clean, prot_clean];
end

%% 4. Run the Triage Diagnostics
report = profile_data(t, X);

%% 5. Plot the Conditioned Trajectories
figure('Name', 'Hes1 Processed Trajectories');
subplot(2,1,1);
plot(t_mrna_raw, mrna_raw, 'k.', 'MarkerSize', 8); hold on;
plot(t, X(:,1), 'b-', 'LineWidth', 2);
ylabel('mRNA Concentration');
legend('Raw Replicates', 'Averaged Trajectory');
grid on;
title('Hes1 Dynamics Diagnostic Visualization');

subplot(2,1,2);
plot(t_prot_raw, prot_raw, 'k.', 'MarkerSize', 8); hold on;
plot(t, X(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Protein Concentration');
legend('Raw Replicates', 'Averaged Trajectory');
grid on;