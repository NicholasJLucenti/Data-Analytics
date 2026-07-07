function samp = check_sampling(t)
%CHECK_SAMPLING Measure sampling interval regularity for a time vector.
%
%   samp = CHECK_SAMPLING(t) returns a struct:
%     .n_points    - number of unique, finite timestamps
%     .dt_mean     - mean sampling interval
%     .dt_std      - std of sampling intervals
%     .dt_cv       - coefficient of variation of dt (std/mean); 0 = perfectly regular
%     .dt_min      - smallest interval (flags near-duplicate samples)
%     .dt_max      - largest interval (flags gaps)
%     .t_span      - total time span (max(t) - min(t))
%     .is_regular  - true if dt_cv is below a small tolerance (default 0.05)
%
%   dt_cv is the key routing signal: near-zero means the data is already
%   on a uniform grid (safe for finite-difference derivatives), while a
%   large value means irregular sampling that favors either densification
%   or a weak/integral-form approach that doesn't require a derivative
%   estimate at all.

t = t(:);
t = t(isfinite(t));
t = unique(sort(t));

samp.n_points = numel(t);

if numel(t) < 2
    samp.dt_mean = NaN;
    samp.dt_std = NaN;
    samp.dt_cv = NaN;
    samp.dt_min = NaN;
    samp.dt_max = NaN;
    samp.t_span = 0;
    samp.is_regular = false;
    return
end

dt = diff(t);
samp.dt_mean = mean(dt);
samp.dt_std = std(dt);
if samp.dt_mean > 0
    samp.dt_cv = samp.dt_std / samp.dt_mean;
else
    samp.dt_cv = NaN;
end
samp.dt_min = min(dt);
samp.dt_max = max(dt);
samp.t_span = t(end) - t(1);

regularityTolerance = 0.05;
samp.is_regular = isfinite(samp.dt_cv) && samp.dt_cv < regularityTolerance;

end