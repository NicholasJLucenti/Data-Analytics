function diagnostics = profile_data(raw)
%PROFILE_DATA Run the full diagnostic suite over a raw data struct array
%and package the results for downstream routing decisions.
%
%   diagnostics = PROFILE_DATA(raw) returns a struct:
%     .per_variable  - struct array, one entry per variable in raw:
%                        .name
%                        .noise     (output of estimate_noise: snr_db, ...)
%                        .sampling  (output of check_sampling: dt_cv, ...)
%                        .sparsity  (output of estimate_sparsity: points_per_period, ...)
%     .horizon       - output of detect_horizon_mismatch across all variables
%     .summary       - pipeline-level rollups consumed by selection/select_route.m:
%                        .min_snr_db        - worst-case SNR across variables (dB)
%                        .max_dt_cv         - worst-case sampling irregularity
%                        .any_sparse        - true if any variable is under-sampled
%                        .overlap_fraction  - shared horizon fraction (from horizon)
%
%   This function assumes raw has already passed io/validate_input.m --
%   it does not re-check structural validity, only computes statistics.

nVars = numel(raw);
if nVars == 0
    error('profile_data:empty', 'raw struct array is empty.');
end

perVar = struct('name', {}, 'noise', {}, 'sampling', {}, 'sparsity', {});

for i = 1:nVars
    if isfield(raw(i), 'name') && ~isempty(raw(i).name)
        nm = raw(i).name;
    else
        nm = sprintf('var%d', i);
    end

    perVar(i).name = nm;
    perVar(i).noise = estimate_noise(raw(i).t, raw(i).y);
    perVar(i).sampling = check_sampling(raw(i).t);
    perVar(i).sparsity = estimate_sparsity(raw(i).t, raw(i).y);
end

diagnostics.per_variable = perVar;
diagnostics.horizon = detect_horizon_mismatch(raw);

snrVals = arrayfun(@(v) v.noise.snr_db, perVar);
cvVals = arrayfun(@(v) v.sampling.dt_cv, perVar);
sparseFlags = arrayfun(@(v) v.sparsity.is_sparse, perVar);

diagnostics.summary.min_snr_db = min(snrVals, [], 'omitnan');
diagnostics.summary.max_dt_cv = max(cvVals, [], 'omitnan');
diagnostics.summary.any_sparse = any(sparseFlags);
diagnostics.summary.overlap_fraction = diagnostics.horizon.overlap_fraction;

end