%% Test Script: Diagnostics on Local Mat-File Hes1 Data
clear; clc; close all;

% 1. Add repository folders to MATLAB path so everything can communicate
addpath('../diagnostics');
addpath('../preprocessing');

fprintf('Loading Hes1 .mat datasets...\n');

%% 2. Load the Data Dynamically
% This safely grabs the first variable out of each binary file
mRNA_struct   = load('../data/Hes1mRNAData.mat');
mRNA_t_struct = load('../data/Hes1mRNATime.mat');
prot_struct   = load('../data/Hes1Data.mat');
prot_t_struct = load('../data/Hes1Time.mat');

mrna_fields = fields(mRNA_struct);
mrna_raw    = mRNA_struct.(mrna_fields{1});

mrna_t_fields = fields(mRNA_t_struct);
t_mrna_raw    = mRNA_t_struct.(mrna_t_fields{1});

prot_fields = fields(prot_struct);
prot_raw    = prot_struct.(prot_fields{1});

prot_t_fields = fields(prot_t_struct);
t_prot_raw    = prot_t_struct.(prot_t_fields{1});

%% 3. Handle Replicates (Averaging multiple measurements at identical time points)
fprintf('Conditioning data: Averaging concurrent replicate measurements...\n');

[t_mrna_uni, ~, idx_mrna] = unique(t_mrna_raw);
mrna_clean = accumarray(idx_mrna, mrna_raw, [], @mean);

[t_prot_uni, ~, idx_prot] = unique(t_prot_raw);
prot_clean = accumarray(idx_prot, prot_raw, [], @mean);

%% 4. Automated Adaptive Alignment & Interpolation
if ~isequal(t_mrna_uni, t_prot_uni)
    fprintf('[INFO] Time grids mismatch. Running automated interpolation selector...\n');
    
    % Call the new function from the preprocessing folder to scan the protein dynamics
    [prot_method, prot_dynamics] = select_interpolation(t_prot_uni, prot_clean);
    fprintf('  -> Protein Dynamics: %s\n', prot_dynamics.reason);
    fprintf('  -> Selected Interpolation Algorithm: "%s"\n', prot_method);
    
    % Call the new function to inspect mRNA dynamics for context
    [mrna_method, mrna_dynamics] = select_interpolation(t_mrna_uni, mrna_clean);
    fprintf('  -> mRNA Dynamics: %s\n', mrna_dynamics.reason);
    
    % Execute the chosen safe interpolation method ('pchip' will be chosen for your protein data)
    prot_aligned = interp1(t_prot_uni, prot_clean, t_mrna_uni, prot_method, 'extrap');
    
    t = t_mrna_uni;
    X = [mrna_clean, prot_aligned];
else
    t = t_mrna_uni;
    X = [mrna_clean, prot_clean];
end

%% 5. Run the Triage Diagnostics Engine
report = profile_data(t, X);

%% 6. Plot the Cleaned and Aligned Trajectories
figure('Name', 'Hes1 Processed Trajectories', 'Position', [100, 100, 800, 500]);

% Top subplot: mRNA
subplot(2,1,1);
plot(t_mrna_raw, mrna_raw, 'k.', 'MarkerSize', 10); hold on;
plot(t, X(:,1), 'b-', 'LineWidth', 2);
ylabel('mRNA Concentration');
legend('Raw Replicates', 'Averaged Trajectory', 'Location', 'best');
grid on;
title('Hes1 Dynamics Adaptive Diagnostic Pipeline');

% Bottom subplot: Protein
subplot(2,1,2);
plot(t_prot_raw, prot_raw, 'k.', 'MarkerSize', 10); hold on;
plot(t, X(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Protein Concentration');
legend('Raw Replicates', 'Averaged Trajectory', 'Location', 'best');
grid on;