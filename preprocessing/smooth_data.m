function X_smooth = smooth_data(t, X, smoothing_factor)
    % SMOOTH_DATA Smooths multi-channel state matrices to suppress noise
    % using a pure MATLAB implementation of Lowess (no toolboxes required).
    %
    % Inputs:
    %   t                - Time vector (N x 1)
    %   X                - Raw state matrix (N x D)
    %   smoothing_factor - Percentage of data points to use in span (e.g., 0.1 to 0.5)

    if nargin < 3
        smoothing_factor = 0.2; 
    end

    [N, D] = size(X);
    X_smooth = zeros(N, D);
    
    % Determine the number of neighbors to consider based on the span factor
    span_pts = max(3, round(smoothing_factor * N));

    for j = 1:D
        xj = X(:, j);
        for i = 1:N
            % 1. Calculate distances from the current point t(i) to all other points
            dists = abs(t - t(i));
            
            % 2. Find the closest k-neighbors
            [sorted_dists, idx] = sort(dists);
            neighbors_idx = idx(1:span_pts);
            max_dist = sorted_dists(span_pts);
            
            % Avoid division by zero if all points are at the same distance
            if max_dist == 0, max_dist = 1; end
            
            % 3. Compute Tricube weights: w = (1 - (d/d_max)^3)^3
            w = (1 - (dists(neighbors_idx) / max_dist).^3).^3;
            w(dists(neighbors_idx) >= max_dist) = 0; % Guardrail boundaries
            
            % 4. Perform locally weighted linear regression
            % Solve: [1, t_local] * beta = x_local with weights W
            W = diag(w);
            A = [ones(span_pts, 1), t(neighbors_idx)];
            y_local = xj(neighbors_idx);
            
            % Weighted least squares formula: beta = (A'*W*A) \ (A'*W*y)
            beta = (A' * W * A) \ (A' * W * y_local);
            
            % 5. Evaluate the local fit at the specific target point t(i)
            X_smooth(i, j) = beta(1) + beta(2) * t(i);
        end
    end
end