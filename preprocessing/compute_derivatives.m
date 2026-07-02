function dXdt = compute_derivatives(t, X)
    % COMPUTE_DERIVATIVES Computes robust time derivatives on regular grids
    % without requiring any MATLAB toolboxes.
    
    [N, D] = size(X);
    dXdt = zeros(N, D);

    for i = 1:D
        % Pure MATLAB native gradient operator (high accuracy on uniform grids)
        dXdt(:, i) = gradient(X(:, i), t);
    end
end