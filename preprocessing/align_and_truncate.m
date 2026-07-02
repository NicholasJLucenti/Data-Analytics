function [t_out, X_out] = align_and_truncate(t_mrna_raw, mrna_raw, t_prot_raw, prot_raw, num_dense_points)
    % ALIGN_AND_TRUNCATE Manages horizon matching and dense uniform resampling.
    %
    % Inputs:
    %   t_mrna_raw       - Raw mRNA time stamps
    %   mrna_raw         - Raw mRNA concentrations
    %   t_prot_raw       - Raw protein time stamps
    %   prot_raw         - Raw protein concentrations
    %   num_dense_points - Target number of densified points (e.g., 300)

    if nargin < 5
        num_dense_points = 300; % Default to healthy density if not specified
    end

    %% 1. Collapse Replicates (Compute means for identical time points)
    [t_mrna_uni, ~, idx_mrna] = unique(t_mrna_raw);
    mrna_clean = accumarray(idx_mrna, mrna_raw, [], @mean);

    [t_prot_uni, ~, idx_prot] = unique(t_prot_raw);
    prot_clean = accumarray(idx_prot, prot_raw, [], @mean);

    %% 2. Horizon Truncation Guardrail
    t_max_common = min(max(t_mrna_uni), max(t_prot_uni));
    t_min_common = max(min(t_mrna_uni), min(t_prot_uni));
    
    % Crop original data vectors to the shared domain boundaries
    idx_mrna_keep = (t_mrna_uni >= t_min_common) & (t_mrna_uni <= t_max_common);
    t_mrna_uni = t_mrna_uni(idx_mrna_keep);
    mrna_clean = mrna_clean(idx_mrna_keep);
    
    idx_prot_keep = (t_prot_uni >= t_min_common) & (t_prot_uni <= t_max_common);
    t_prot_uni = t_prot_uni(idx_prot_keep);
    prot_clean = prot_clean(idx_prot_keep);

    %% 3. Create the Dense, Uniform Time Grid
    % This generates 300 perfectly spaced time points between the boundaries
    t_out = linspace(t_min_common, t_max_common, num_dense_points)';
    fprintf('[DENSIFICATION] Resampling %d sparse points into %d uniform steps.\n', ...
            length(t_mrna_uni), num_dense_points);

    %% 4. Adaptive Interpolation onto the Dense Grid
    % Profile and interpolate mRNA
    [mrna_method, mrna_dynamics] = select_interpolation(t_mrna_uni, mrna_clean, t_out);
    if strcmp(mrna_method, 'fourier')
        fit_mrna = fit(t_mrna_uni, mrna_clean, 'fourier1');
        mrna_dense = fit_mrna(t_out);
    else
        mrna_dense = interp1(t_mrna_uni, mrna_clean, t_out, mrna_method, 'extrap');
    end

    % Profile and interpolate Protein
    [prot_method, prot_dynamics] = select_interpolation(t_prot_uni, prot_clean, t_out);
    if strcmp(prot_method, 'fourier')
        fit_prot = fit(t_prot_uni, prot_clean, 'fourier1');
        prot_dense = fit_prot(t_out);
    else
        prot_dense = interp1(t_prot_uni, prot_clean, t_out, prot_method, 'extrap');
    end

    X_out = [mrna_dense, prot_dense];
end