function [Xi, library_names, L, Y, weak_info] = run_weak_sindy(t, X, lambda, poly_order, varargin)
%RUN_WEAK_SINDY Discover governing equations via weak/integral-form SINDy.
%
%   [Xi, library_names, L, Y, weak_info] = RUN_WEAK_SINDY(t, X, lambda, poly_order, ...)
%
%   Unlike standard SINDy, this variant never estimates a pointwise
%   derivative of the (possibly noisy) raw signal. Instead it multiplies
%   the ODE dx/dt = Theta(x)*xi by a smooth, compactly-supported test
%   function phi and integrates by parts over a local window:
%
%       d/dt x(t) = Theta(x(t)) * xi
%       ... multiply by phi, integrate over [t_k-r, t_k+r], integrate
%           the left side by parts (phi vanishes at the window edges,
%           so the boundary term drops out) ...
%       -integral( phi'(t) x(t) dt ) = integral( phi(t) Theta(x(t)) dt ) * xi
%
%   Each test-function window k contributes one row to a linear system
%   A*xi = b, built entirely from integrals of the (untouched) signal and
%   the library evaluated on it -- no differentiation of noisy data
%   anywhere. This is what makes weak-form SINDy robust on exactly the
%   kind of sparse, low-SNR data that breaks standard SINDy's finite
%   differences (see selection/select_route.m's SNR/sparsity gates).
%
%   Inputs:
%     t          - dense, uniformly-spaced time vector (N x 1). Use
%                   preprocessing/align_and_truncate.m's output directly --
%                   do NOT smooth X before calling this; weak-form
%                   integration is the noise-handling step.
%     X          - state matrix (N x D)
%     lambda     - sparsity threshold (default 0.1). NOTE: the natural
%                   scale of the weak-form linear system (A, b) differs
%                   from the pointwise (Theta, dXdt) system, since entries
%                   are integrals over a window rather than instantaneous
%                   values -- lambda for this variant will typically need
%                   separate tuning from the standard-variant lambda.
%     poly_order - highest polynomial degree in the candidate library
%                   (default 2)
%
%   Name-value options:
%     'NumTestFunctions' - target number of test-function windows (default 150)
%     'WindowPoints'     - number of grid points spanned by each test
%                           function's support; forced odd (default 21)
%     'TestFunctionOrder' - smoothness exponent p of the bump test
%                           function phi(s) = (1-s^2)^p (default 4)
%
%   Outputs:
%     Xi            - sparse coefficient matrix (library_size x D)
%     library_names - cell array describing what each row represents
%     L             - the weak-form library matrix actually solved
%                     against (= A, one row per test function). Returned
%                     so selection/ensemble_sindy.m can bootstrap over
%                     test functions without rebuilding them.
%     Y             - the weak-form target matrix actually solved against
%                     (= b). Same rationale as L.
%     weak_info     - struct with num_test_functions, window_points,
%                     radius, test_function_order, for diagnostics/reporting

p = inputParser;
addParameter(p, 'NumTestFunctions', 150, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'WindowPoints', 21, @(x) isnumeric(x) && isscalar(x) && x >= 5);
addParameter(p, 'TestFunctionOrder', 4, @(x) isnumeric(x) && isscalar(x) && x >= 2);
parse(p, varargin{:});
opts = p.Results;

if nargin < 4 || isempty(poly_order)
    poly_order = 2;
end
if nargin < 3 || isempty(lambda)
    lambda = 0.1;
end

t = t(:);
[N, D] = size(X);
if numel(t) ~= N
    error('run_weak_sindy:sizeMismatch', 'numel(t) must equal size(X,1).');
end
if N < 10
    error('run_weak_sindy:tooFewPoints', 'Need at least 10 time points to place test functions.');
end

dt = median(diff(t));
if dt <= 0
    error('run_weak_sindy:badTimeVector', 't must be strictly increasing.');
end
if std(diff(t)) / dt > 1e-3
    warning('run_weak_sindy:nonUniformGrid', ...
        ['Time grid does not look uniform; weak-form quadrature here assumes uniform spacing. ' ...
         'Run preprocessing/align_and_truncate.m first to densify onto a uniform grid.']);
end

%% 1. Set up test-function windows
windowPoints = min(opts.WindowPoints, 2*floor((N-1)/2) - 1); % keep inside domain
windowPoints = max(windowPoints, 5);
if mod(windowPoints, 2) == 0
    windowPoints = windowPoints + 1; % force odd for a symmetric window
end
halfWin = (windowPoints - 1) / 2;
r = halfWin * dt; % physical radius of each test function's support

minIdx = halfWin + 1;
maxIdx = N - halfWin;
if maxIdx <= minIdx
    error('run_weak_sindy:domainTooSmall', ...
        'Not enough points to place any test functions with WindowPoints=%d (N=%d).', windowPoints, N);
end

numTF = min(opts.NumTestFunctions, maxIdx - minIdx + 1);
centerIdx = unique(round(linspace(minIdx, maxIdx, numTF)));
numTF = numel(centerIdx);

%% 2. Build the full candidate library once, evaluated at every grid point
[Theta_full, library_names] = build_library(X, poly_order);
M = size(Theta_full, 2);

%% 3. Assemble the weak-form linear system: one row per test function
A = zeros(numTF, M);
b = zeros(numTF, D);
p_order = opts.TestFunctionOrder;

for k = 1:numTF
    ci = centerIdx(k);
    idxRange = (ci - halfWin):(ci + halfWin);
    tWin = t(idxRange);
    tk = t(ci);
    s = (tWin - tk) / r;

    phi = (1 - s.^2).^p_order;
    phi(abs(s) >= 1) = 0;
    dphi = (-2 * p_order / r) .* s .* (1 - s.^2).^(p_order - 1);
    dphi(abs(s) >= 1) = 0;

    ThetaWin = Theta_full(idxRange, :);
    XWin = X(idxRange, :);

    for m = 1:M
        A(k, m) = trapz(tWin, phi .* ThetaWin(:, m));
    end
    for d = 1:D
        b(k, d) = -trapz(tWin, dphi .* XWin(:, d));
    end
end

L = A;
Y = b;

%% 4. Sparse regression on the integrated system
Xi = stlsq_solve(L, Y, lambda, 10);

weak_info.num_test_functions = numTF;
weak_info.window_points = windowPoints;
weak_info.radius = r;
weak_info.test_function_order = p_order;

end