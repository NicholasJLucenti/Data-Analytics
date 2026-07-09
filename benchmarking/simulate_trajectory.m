function [t_sim, X_sim, sim_info] = simulate_trajectory(model, poly_order, t_span, x0, t_eval, varargin)
%SIMULATE_TRAJECTORY Forward-integrate a discovered SINDy model, of any variant.
%
%   [t_sim, X_sim, sim_info] = SIMULATE_TRAJECTORY(model, poly_order, t_span, x0, t_eval, ...)
%
%   Inputs:
%     model      - EITHER a sparse coefficient matrix (M x D) from
%                  run_standard_sindy.m/run_weak_sindy.m (dxdt = Theta(x)*model),
%                  OR a 1xD struct array from run_implicit_sindy.m with
%                  .numerator_Xi/.denominator_Xi fields (dxdt(d) =
%                  Theta(x)*numerator_Xi(d) / Theta(x)*denominator_Xi(d)).
%                  Dispatch is automatic based on whether this is numeric
%                  or a struct.
%     poly_order - polynomial order used to build the library model was fit
%                  against (must match, since the RHS is reconstructed
%                  from libraries/build_polynomial_library.m at every step)
%     t_span     - [t_start, t_end]
%     x0         - initial condition (1 x D or D x 1)
%     t_eval     - (optional) vector of times at which to report the
%                  solution, so it lines up directly with real data for
%                  comparison. If omitted, ode45 chooses its own steps.
%
%   Name-value options:
%     'DivergenceBound' - if any state exceeds this magnitude, integration
%                          stops immediately via a terminal ODE event,
%                          instead of letting the adaptive solver grind
%                          through ever-smaller steps trying to resolve a
%                          near-singular blow-up to the requested
%                          tolerance (default 1e6, matches the divergence
%                          check below)
%     'MaxWallSeconds'   - hard wall-clock cap on a single simulation call
%                          (default 3). This is a backstop for cases that
%                          aren't diverging past DivergenceBound but are
%                          just stiff/slow for a bounded, non-diverging
%                          trajectory -- without this, one pathological
%                          candidate in a large hyperparameter sweep can
%                          stall the whole run.
%
%   Outputs:
%     t_sim    - times of the returned solution (== t_eval if given; will
%                be SHORTER than requested if a divergence event or the
%                watchdog cut integration off early)
%     X_sim    - simulated trajectory (numel(t_sim) x D)
%     sim_info - struct: .success (bool), .message (string, e.g. why it
%                failed, diverged, or was cut off)
%
%   A model that diverges or blows up under forward simulation is exactly
%   the failure mode that in-sample regression fit (AIC/BIC, R^2, even
%   ensemble inclusion probability) cannot see -- this is the check that
%   catches it. Cutting such cases off early (rather than fighting to
%   integrate through them precisely) is correct here: we only need to
%   know that the model diverged, not exactly how.

p = inputParser;
addParameter(p, 'DivergenceBound', 1e6, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MaxWallSeconds', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});
opts = p.Results;

x0 = x0(:)';
sim_info.success = true;
sim_info.message = '';

startTime = tic;
rhs = @(tt, xx) local_rhs(xx, model, poly_order, startTime, opts.MaxWallSeconds);
odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
    'Events', @(tt, xx) local_divergence_event(xx, opts.DivergenceBound), ...
    'OutputFcn', @(tt, xx, flag) local_watchdog(tt, xx, flag, startTime, opts.MaxWallSeconds));

% Unstable candidate models (common across a wide hyperparameter sweep)
% cause ode45 to report a "failed to meet tolerances" console warning
% when it hits the terminal event / watchdog cutoff. That failure is
% already detected below via truncated/non-finite output, so the warning
% itself is just noise across a large sweep -- suppress it here,
% restoring the previous warning state automatically even if this
% function errors out.
warnState = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>

try
    if nargin >= 5 && ~isempty(t_eval)
        [t_sim, X_sim] = ode45(rhs, t_eval, x0, odeOpts);
    else
        [t_sim, X_sim] = ode45(rhs, t_span, x0, odeOpts);
    end

    if isempty(X_sim) || any(~isfinite(X_sim(:))) || any(abs(X_sim(:)) > opts.DivergenceBound)
        sim_info.success = false;
        sim_info.message = 'Trajectory diverged (non-finite or unbounded values encountered).';
    elseif nargin >= 5 && ~isempty(t_eval) && numel(t_sim) < numel(t_eval)
        sim_info.success = false;
        if toc(startTime) >= opts.MaxWallSeconds
            sim_info.message = sprintf('Integration cut off by wall-clock watchdog (> %.1fs).', opts.MaxWallSeconds);
        else
            sim_info.message = 'Integration terminated early (divergence event triggered before reaching requested end time).';
        end
    end
catch ME
    t_sim = [];
    X_sim = [];
    sim_info.success = false;
    sim_info.message = sprintf('Integration failed: %s', ME.message);
end

end


function dxdt = local_rhs(x, model, poly_order, startTime, maxSeconds)
    % model is either:
    %   - numeric (M x D) matrix: polynomial variant (standard/weak) --
    %     dxdt = Theta(x)*model
    %   - struct array (1 x D) with .numerator_Xi/.denominator_Xi: implicit
    %     variant -- dxdt(d) = Theta(x)*numerator_Xi(d) / Theta(x)*denominator_Xi(d)
    %
    % This is called on EVERY derivative evaluation ode45 makes, including
    % internal step-size retries that never become an accepted step (and
    % therefore never trigger OutputFcn/local_watchdog below). A stiff or
    % near-singular RHS -- e.g. a Hill denominator crossing zero when an
    % unstable trajectory drives a state negative with an odd exponent --
    % can make ode45 spend a long time re-attempting a single step at
    % ever-finer resolution without ever calling OutputFcn again. Checking
    % the wall-clock budget and rejecting non-finite output HERE, not just
    % in OutputFcn, is what actually stops that stall: throwing aborts
    % ode45 immediately rather than letting it grind.
    if toc(startTime) > maxSeconds
        error('simulate_trajectory:rhsWatchdogTimeout', ...
            'RHS wall-clock budget (%.1fs) exceeded mid-step (likely stuck retrying near a singularity).', maxSeconds);
    end

    xrow = x(:)'; % build_library expects (N x D); N=1 for a single state
    Theta_row = build_library(xrow, poly_order);

    if isstruct(model)
        D = numel(model);
        dxdt = zeros(D, 1);
        for d = 1:D
            num_val = Theta_row * model(d).numerator_Xi;
            den_val = Theta_row * model(d).denominator_Xi;
            dxdt(d) = num_val / den_val;
        end
    else
        dxdt = (Theta_row * model)';
    end

    if ~all(isfinite(dxdt))
        error('simulate_trajectory:nonFiniteDerivative', ...
            'Non-finite derivative encountered (likely a near-zero denominator or unbounded state).');
    end
end


function [value, isterminal, direction] = local_divergence_event(x, bound)
    % Fires (and stops integration) the moment any state variable's
    % magnitude crosses the divergence bound, rather than letting ode45
    % try to precisely resolve a trajectory that's headed to infinity.
    value = bound - max(abs(x));
    isterminal = 1;
    direction = -1;
end


function status = local_watchdog(~, ~, flag, startTime, maxSeconds)
    % ode45 OutputFcn: returning status=1 halts integration immediately.
    % Called with flag='init'/'done' at the start/end and '' at every
    % accepted step in between -- only check the clock on the per-step
    % calls.
    status = 0;
    if isempty(flag) && toc(startTime) > maxSeconds
        status = 1;
    end
end