function raw = load_raw_data(source, varargin)
%LOAD_RAW_DATA Ingest raw experimental data into the pipeline's standard format.
%
%   raw = LOAD_RAW_DATA(source) accepts:
%     - a path to a .mat file containing one struct-per-variable, or a
%       single struct with fields like t_mRNA, y_mRNA, t_protein, y_protein
%     - a path to a folder containing multiple .mat files (one per channel)
%     - an in-memory struct/struct array already in raw form (fields t, y)
%
%   Returns:
%     raw - 1xN struct array, one entry per state variable, with fields:
%       .name   - string identifier (e.g. 'mRNA')
%       .t      - column vector of raw timestamps (may contain duplicates/replicates)
%       .y      - column vector of raw observations, same length as t
%       .source - string describing where this variable came from (provenance)
%
%   Name-value options:
%     'VariableNames' - cell array of names to assign in order, overrides
%                        whatever names were auto-detected.
%
%   This function does NOT clean, truncate, or validate the data -- it
%   only normalizes structure so every downstream stage can rely on a
%   consistent shape. See validate_input.m for sanity checks and
%   preprocessing/align_and_truncate.m for cleaning.

p = inputParser;
addParameter(p, 'VariableNames', {}, @iscell);
parse(p, varargin{:});
varNames = p.Results.VariableNames;

if isstruct(source)
    raw = local_parse_struct(source, varNames);
elseif ischar(source) || isstring(source)
    source = char(source);
    if isfolder(source)
        raw = local_parse_folder(source, varNames);
    elseif isfile(source)
        raw = local_parse_matfile(source, varNames);
    else
        error('load_raw_data:invalidSource', ...
            'Source path "%s" is neither a valid file nor folder.', source);
    end
else
    error('load_raw_data:invalidSource', ...
        'Source must be a struct, a .mat file path, or a folder path.');
end

if isempty(raw)
    error('load_raw_data:empty', 'No variables could be parsed from the provided source.');
end

end


% ---------------------------------------------------------------------
% local helpers
% ---------------------------------------------------------------------

function raw = local_parse_matfile(filepath, varNames)
    s = load(filepath);
    raw = local_parse_struct(s, varNames);
    for i = 1:numel(raw)
        raw(i).source = filepath;
    end
end

function raw = local_parse_folder(folderpath, varNames)
    files = dir(fullfile(folderpath, '*.mat'));
    if isempty(files)
        error('load_raw_data:noMatFiles', 'No .mat files found in folder "%s".', folderpath);
    end
    raw = struct('name', {}, 't', {}, 'y', {}, 'source', {});
    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);
        s = load(fpath);
        entry = local_parse_struct(s, {});
        for j = 1:numel(entry)
            entry(j).source = fpath;
            if isempty(entry(j).name) || strcmp(entry(j).name, sprintf('var%d', j))
                [~, base, ~] = fileparts(files(i).name);
                entry(j).name = base;
            end
        end
        raw = [raw, entry]; %#ok<AGROW>
    end
    if ~isempty(varNames)
        for i = 1:min(numel(varNames), numel(raw))
            raw(i).name = varNames{i};
        end
    end
end

function raw = local_parse_struct(s, varNames)
    % Case 1: s already looks like a raw struct array with t/y fields
    if isfield(s, 't') && isfield(s, 'y')
        raw = s;
        for i = 1:numel(raw)
            raw(i).t = double(raw(i).t(:));
            raw(i).y = double(raw(i).y(:));
            if ~isfield(raw(i), 'name') || isempty(raw(i).name)
                raw(i).name = sprintf('var%d', i);
            end
            if ~isfield(raw(i), 'source')
                raw(i).source = '';
            end
        end
        if ~isempty(varNames)
            for i = 1:min(numel(varNames), numel(raw))
                raw(i).name = varNames{i};
            end
        end
        return
    end

    % Case 2: s is a container (e.g. loaded .mat) with fields like
    % t_mRNA, y_mRNA, t_protein, y_protein  -> auto-pair by suffix
    fn = fieldnames(s);
    tFields = fn(startsWith(fn, 't_') | strcmpi(fn, 't'));
    raw = struct('name', {}, 't', {}, 'y', {}, 'source', {});
    count = 0;
    for i = 1:numel(tFields)
        tf = tFields{i};
        if strcmpi(tf, 't')
            suffix = '';
        else
            suffix = extractAfter(tf, 't_');
        end
        if isempty(suffix)
            yf = 'y';
        else
            yf = ['y_' suffix];
        end
        if isfield(s, yf)
            count = count + 1;
            raw(count).name = suffix;
            if isempty(raw(count).name)
                raw(count).name = sprintf('var%d', count);
            end
            raw(count).t = double(s.(tf)(:));
            raw(count).y = double(s.(yf)(:));
            raw(count).source = '';
        end
    end

    if count == 0
        error('load_raw_data:unrecognizedFormat', ...
            ['Could not auto-detect variable pairs. Expected fields named ' ...
             '"t"/"y" or "t_<name>"/"y_<name>", or use the VariableNames option.']);
    end

    if ~isempty(varNames)
        for i = 1:min(numel(varNames), numel(raw))
            raw(i).name = varNames{i};
        end
    end
end