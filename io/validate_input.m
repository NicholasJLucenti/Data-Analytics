function report = validate_input(raw)
%VALIDATE_INPUT Sanity-check a raw data struct array before it enters the pipeline.
%
%   report = VALIDATE_INPUT(raw) runs structural and numerical checks on
%   each variable in raw (as produced by load_raw_data) and returns:
%     .ok        - true if no critical issues were found
%     .variables - 1xN struct array, one per raw variable:
%                    .name, .n_points, .issues (cell array), .warnings (cell array)
%     .messages  - flat cell array of all issue/warning strings, for quick display
%
%   Critical issues (report.ok = false) include: missing t/y fields,
%   mismatched t/y lengths, fewer than 3 finite points. These should stop
%   the pipeline. Non-critical problems (duplicate timestamps, unsorted
%   time, isolated NaNs, zero variance) are recorded as warnings and left
%   for downstream stages (which already handle replicates/truncation)
%   to resolve -- they are surfaced here so the user knows they happened.

nVars = numel(raw);
report.ok = true;
report.variables = struct('name', {}, 'n_points', {}, 'issues', {}, 'warnings', {});
report.messages = {};

if nVars == 0
    report.ok = false;
    report.messages{end+1} = 'no variables found in raw input';
    return
end

for i = 1:nVars
    v = raw(i);
    issues = {};
    warnings = {};
    nm = local_name(v, i);

    if ~isfield(v, 't') || ~isfield(v, 'y') || isempty(v.t) || isempty(v.y)
        issues{end+1} = 'missing or empty t/y field'; %#ok<AGROW>
        report.variables(i) = struct('name', nm, 'n_points', 0, ...
            'issues', {issues}, 'warnings', {warnings});
        report.ok = false;
        report.messages{end+1} = sprintf('[%s] ERROR: %s', nm, issues{1});
        continue
    end

    t = v.t(:); y = v.y(:);

    if numel(t) ~= numel(y)
        issues{end+1} = sprintf('length mismatch: numel(t)=%d, numel(y)=%d', numel(t), numel(y)); %#ok<AGROW>
    end

    finiteMask = isfinite(t) & isfinite(y);
    nFinite = sum(finiteMask);
    if nFinite < 3
        issues{end+1} = sprintf('fewer than 3 finite (t,y) pairs (found %d)', nFinite); %#ok<AGROW>
    end

    nBadY = sum(~isfinite(y));
    if nBadY > 0 && nBadY < numel(y)
        warnings{end+1} = sprintf('%d non-finite y values will be dropped downstream', nBadY); %#ok<AGROW>
    end

    if numel(t) > 1
        nDup = numel(t) - numel(unique(t));
        if nDup > 0
            warnings{end+1} = sprintf('%d duplicate timestamps detected (will be treated as replicates)', nDup); %#ok<AGROW>
        end
        if ~issorted(t)
            warnings{end+1} = 'timestamps are not sorted ascending'; %#ok<AGROW>
        end
    end

    if nFinite >= 2 && var(y(finiteMask)) == 0
        warnings{end+1} = 'signal has zero variance (constant) over all finite points'; %#ok<AGROW>
    end

    report.variables(i) = struct('name', nm, 'n_points', numel(t), ...
        'issues', {issues}, 'warnings', {warnings});

    if ~isempty(issues)
        report.ok = false;
    end

    for k = 1:numel(issues)
        report.messages{end+1} = sprintf('[%s] ERROR: %s', nm, issues{k}); %#ok<AGROW>
    end
    for k = 1:numel(warnings)
        report.messages{end+1} = sprintf('[%s] warning: %s', nm, warnings{k}); %#ok<AGROW>
    end
end

end


function n = local_name(v, i)
    if isfield(v, 'name') && ~isempty(v.name)
        n = v.name;
    else
        n = sprintf('var%d', i);
    end
end