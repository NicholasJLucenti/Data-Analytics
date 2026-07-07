function [t_out, X_out, names] = align_and_truncate(raw, num_dense_points)
%ALIGN_AND_TRUNCATE Collapse replicates, truncate to the shared time
%horizon, and densify an arbitrary number of raw channels onto one
%uniform grid. This is the seam between the pre-merge struct-array
%convention (io/, diagnostics/) and the merged (t, X) matrix convention
%used by everything downstream (libraries/, variants/, selection/).
%
%   [t_out, X_out, names] = ALIGN_AND_TRUNCATE(raw, num_dense_points)
%
%   Inputs:
%     raw               - struct array from io/load_raw_data.m, with
%                          fields .name, .t, .y (any number of channels)
%     num_dense_points  - target number of uniform output points
%                          (default 300)
%
%   Outputs:
%     t_out  - shared dense time vector (num_dense_points x 1)
%     X_out  - densified state matrix (num_dense_points x D), column
%              order matches raw
%     names  - 1xD cell array of channel names, matching X_out columns

if nargin < 2 || isempty(num_dense_points)
    num_dense_points = 300;
end

nVars = numel(raw);
if nVars == 0
    error('align_and_truncate:empty', 'raw struct array is empty.');
end

%% 1. Collapse replicates (mean of identical timestamps) per channel
tClean = cell(nVars, 1);
yClean = cell(nVars, 1);
for i = 1:nVars
    t = raw(i).t(:); y = raw(i).y(:);
    valid = isfinite(t) & isfinite(y);
    t = t(valid); y = y(valid);
    [tu, ~, ic] = unique(t);
    yu = accumarray(ic, y, [], @mean);
    tClean{i} = tu;
    yClean{i} = yu;
end

%% 2. Horizon truncation guardrail (shared across all channels)
rawForHorizon = raw;
for i = 1:nVars
    rawForHorizon(i).t = tClean{i};
    rawForHorizon(i).y = yClean{i};
end
horizon = detect_horizon_mismatch(rawForHorizon);

if horizon.t_end <= horizon.t_start
    error('align_and_truncate:noOverlap', ...
        'Channels do not share a common time window; cannot align.');
end

for i = 1:nVars
    keep = tClean{i} >= horizon.t_start & tClean{i} <= horizon.t_end;
    tClean{i} = tClean{i}(keep);
    yClean{i} = yClean{i}(keep);
    if numel(tClean{i}) < 2
        error('align_and_truncate:tooFewPointsAfterTruncation', ...
            'Channel %d has fewer than 2 points left after horizon truncation.', i);
    end
end

%% 3. Dense uniform target grid
t_out = linspace(horizon.t_start, horizon.t_end, num_dense_points)';
fprintf('[DENSIFICATION] Mapping %d channel(s) onto %d uniform steps over [%.4g, %.4g].\n', ...
    nVars, num_dense_points, horizon.t_start, horizon.t_end);

%% 4. Adaptive interpolation per channel
X_out = zeros(num_dense_points, nVars);
names = cell(1, nVars);
for i = 1:nVars
    if isfield(raw(i), 'name') && ~isempty(raw(i).name)
        names{i} = raw(i).name;
    else
        names{i} = sprintf('var%d', i);
    end

    [method, report] = select_interpolation(tClean{i}, yClean{i}, t_out);
    fprintf('  -> %s: interpolation = %s (%s)\n', names{i}, method, report.reason);

    if strcmp(method, 'fourier')
        fit_fn = fit(tClean{i}, yClean{i}, 'fourier1');
        X_out(:, i) = fit_fn(t_out);
    else
        X_out(:, i) = interp1(tClean{i}, yClean{i}, t_out, method, 'extrap');
    end
end

end