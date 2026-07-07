function [Theta, names] = build_polynomial_library(X, order)
%BUILD_POLYNOMIAL_LIBRARY Generate polynomial candidate terms up to a
%given order for SINDy-style sparse regression.
%
%   [Theta, names] = BUILD_POLYNOMIAL_LIBRARY(X, order)
%
%   Inputs:
%     X     - State matrix (N x D)
%     order - Highest polynomial degree to include (0, 1, 2, or 3).
%             Default 2.
%
%   Outputs:
%     Theta - Candidate library matrix (N x number_of_terms)
%     names - Cell array of human-readable term names, in the same
%             column order as Theta
%
%   Extracted from variants/run_standard_sindy.m so other variants and
%   selection/evaluate_pareto.m can reuse (and eventually augment, e.g.
%   with rational terms) the same library-building logic.

if nargin < 2 || isempty(order)
    order = 2;
end

[N, D] = size(X);

Theta = ones(N, 1);
names = {'1'};

if order >= 1
    Theta = [Theta, X];
    for i = 1:D
        names{end+1} = sprintf('x%d', i); %#ok<AGROW>
    end
end

if order >= 2
    for i = 1:D
        for j = i:D
            Theta = [Theta, X(:,i).*X(:,j)]; %#ok<AGROW>
            names{end+1} = sprintf('x%d*x%d', i, j); %#ok<AGROW>
        end
    end
end

if order >= 3
    for i = 1:D
        for j = i:D
            for k = j:D
                Theta = [Theta, X(:,i).*X(:,j).*X(:,k)]; %#ok<AGROW>
                names{end+1} = sprintf('x%d*x%d*x%d', i, j, k); %#ok<AGROW>
            end
        end
    end
end

end