function print_equations(model, library_names, state_names)
%PRINT_EQUATIONS Print a discovered SINDy model in readable form,
%regardless of which variant produced it.
%
%   PRINT_EQUATIONS(model, library_names, state_names)
%
%   Inputs:
%     model         - coefficient matrix (M x D) from standard/weak
%                     SINDy, OR a 1xD rational struct array (with
%                     .numerator_Xi/.denominator_Xi) from run_implicit_sindy.m
%     library_names - candidate term names matching model's rows
%     state_names   - cell array of state variable names (for the d/dt labels)

if isstruct(model)
    for d = 1:numel(model)
        num_str = local_format_poly(model(d).numerator_Xi, library_names);
        den_str = local_format_poly(model(d).denominator_Xi, library_names);
        fprintf('  d%s/dt = (%s) / (%s)\n', state_names{d}, num_str, den_str);
    end
else
    for d = 1:size(model, 2)
        fprintf('  d%s/dt = %s\n', state_names{d}, local_format_poly(model(:, d), library_names));
    end
end

end


function s = local_format_poly(coeffs, names)
    active = find(coeffs ~= 0);
    if isempty(active)
        s = '0';
        return
    end
    parts = cell(1, numel(active));
    for k = 1:numel(active)
        idx = active(k);
        parts{k} = sprintf('%+.4f*%s', coeffs(idx), names{idx});
    end
    s = strtrim(strjoin(parts, ' '));
end