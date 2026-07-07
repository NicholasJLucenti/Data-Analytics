function horizon = detect_horizon_mismatch(raw)
%DETECT_HORIZON_MISMATCH Find the shared time horizon across multiple
%state-variable channels that may start/stop at different times.
%
%   horizon = DETECT_HORIZON_MISMATCH(raw) takes a raw struct array (as
%   produced by load_raw_data) and returns a struct:
%     .t_start          - latest of the individual start times (max of mins)
%     .t_end            - earliest of the individual end times (min of maxs)
%     .per_variable      - struct array with .name, .t_min, .t_max, .points_dropped
%     .has_mismatch      - true if variables do not all share the same start/end
%     .overlap_fraction  - fraction of the full observed range
%                          [global min, global max] that survives
%                          truncation to [t_start, t_end]
%
%   This is a pure diagnostic -- it does not modify data. It mirrors the
%   "temporal truncation guardrail" logic used in
%   preprocessing/align_and_truncate.m, factored out here so the
%   truncation decision can be inspected/logged before it's applied.

nVars = numel(raw);
if nVars == 0
    error('detect_horizon_mismatch:empty', 'raw struct array is empty.');
end

perVar = struct('name', {}, 't_min', {}, 't_max', {}, 'points_dropped', {});
starts = nan(1, nVars);
ends = nan(1, nVars);

for i = 1:nVars
    t = raw(i).t(:);
    t = t(isfinite(t));
    if isempty(t)
        starts(i) = NaN;
        ends(i) = NaN;
    else
        starts(i) = min(t);
        ends(i) = max(t);
    end
    if isfield(raw(i), 'name') && ~isempty(raw(i).name)
        nm = raw(i).name;
    else
        nm = sprintf('var%d', i);
    end
    perVar(i).name = nm;
    perVar(i).t_min = starts(i);
    perVar(i).t_max = ends(i);
    perVar(i).points_dropped = 0; % filled in below
end

if any(isnan(starts)) || any(isnan(ends))
    warning('detect_horizon_mismatch:emptyChannel', ...
        'One or more variables have no finite timestamps; horizon excludes them.');
end

horizon.t_start = max(starts, [], 'omitnan');
horizon.t_end = min(ends, [], 'omitnan');

globalMin = min(starts, [], 'omitnan');
globalMax = max(ends, [], 'omitnan');

for i = 1:nVars
    t = raw(i).t(:);
    t = t(isfinite(t));
    dropped = sum(t < horizon.t_start | t > horizon.t_end);
    perVar(i).points_dropped = dropped;
end

horizon.per_variable = perVar;
horizon.has_mismatch = ~(all(starts == starts(1)) && all(ends == ends(1)));

fullRange = globalMax - globalMin;
sharedRange = max(horizon.t_end - horizon.t_start, 0);
if fullRange > 0
    horizon.overlap_fraction = sharedRange / fullRange;
else
    horizon.overlap_fraction = 1;
end

if horizon.t_end <= horizon.t_start
    warning('detect_horizon_mismatch:noOverlap', ...
        'Variables do not share a common time overlap; downstream truncation will fail.');
end

end