function [Xi, library_names] = run_standard_sindy(t, X, lambda, poly_order)
    % RUN_STANDARD_SINDY disovers governing equations using classic STLSQ.
    %
    % Inputs:
    %   t          - Truncated time vector (N x 1)
    %   X          - Truncated state matrix (N x D) -> [mRNA, Protein]
    %   lambda     - Sparsity threshold hyperparameter (e.g., 0.1)
    %   poly_order - Highest polynomial degree for library (e.g., 2 or 3)
    %
    % Outputs:
    %   Xi            - Sparse coefficient matrix (Library_Size x D)
    %   library_names - Cell array describing what each row represents

    [N, D] = size(X);

    %% 1. Compute Numerical Derivatives (dX/dt)
  % New way using your modular preprocessing files!
    X_clean = smooth_data(t, X, 0.4); % Step 1: Filter out lingering jitters
    dXdt = compute_derivatives(t, X_clean); % Step 2: Extract robust gradients

    %% 2. Construct Candidate Library Theta(X)
    % For a 2D system up to 2nd order, this creates: [1, x1, x2, x1^2, x1*x2, x2^2]
    [Theta, library_names] = build_polynomial_library(X, poly_order);

    %% 3. Sequential Thresholded Least Squares (STLSQ) Optimizer
    % This is the heart of the standard SINDy paradigm
    Library_Size = size(Theta, 2);
    Xi = Theta \ dXdt; % Initial unthresholded least-squares guess

    % Loop to enforce sparsity
    for iter = 1:10
        small_indices = abs(Xi) < lambda; % Find coefficients smaller than lambda
        Xi(small_indices) = 0;             % Force them to zero
        
        for j = 1:D
            big_indices = ~small_indices(:, j);
            % Re-regress over only the remaining active terms
            Xi(big_indices, j) = Theta(:, big_indices) \ dXdt(:, j);
        end
    end
end

function [Theta, names] = build_polynomial_library(X, order)
    % Helper function to generate polynomial combinations up to 'order'
    [N, D] = size(X);
    
    % Order 0 (Constant term)
    Theta = ones(N, 1);
    names = {'1'};
    
    % Order 1 (Linear terms)
    if order >= 1
        Theta = [Theta, X];
        for i = 1:D
            names{end+1} = sprintf('x%d', i);
        end
    end
    
    % Order 2 (Quadratic interactions)
    if order >= 2
        for i = 1:D
            for j = i:D
                Theta = [Theta, X(:,i).*X(:,j)];
                names{end+1} = sprintf('x%d*x%d', i, j);
            end
        end
    end
    
    % Order 3 (Cubic interactions)
    if order >= 3
        for i = 1:D
            for j = i:D
                for k = j:D
                    Theta = [Theta, X(:,i).*X(:,j).*X(:,k)];
                    names{end+1} = sprintf('x%d*x%d*x%d', i, j, k);
                end
            end
        end
    end
end