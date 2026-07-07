function [t_sim, X_sim, sim_info] = simulate_trajectory(Xi, poly_order, t_span, x0, t_eval)
%SIMULATE_TRAJECTORY Forward-integrate a discovered polynomial ODE model.
%
%   [t_sim, X_sim, sim_info] = SIMULATE_TRAJECTORY(Xi, poly_order, t_span, x0, t_eval)
%
%   Inputs:
%     Xi         - sparse coefficient matrix (M x D) from a SINDy variant
%     poly_order - polynomial order used to build the library Xi was fit
%                  against (must match, since the RHS is reconstructed
%                  from libraries/build_polynomial_library.m at every step)
%     t_span     - [t_start, t_end]
%     x0         - initial condition (1 x D or D x 1)
%     t_eval     - (optional) vector of times at which to report the
%                  solution, so it lines up directly with real data for
%                  comparison. If omitted, ode45 chooses its own steps.
%
%   Outputs:
%     t_sim    - times of the returned solution (== t_eval if given)
%     X_sim    - simulated trajectory (numel(t_sim) x D)
%     sim_info - struct: .success (bool), .message (string, e.g. why it
%                failed or diverged)
%
%   A model that diverges or blows up under forward simulation is exactly
%   the failure mode that in-sample regression fit (AIC/BIC, R^2, even
%   ensemble inclusion probability) cannot see -- this is the check that
%   catches it.

x0 = x0(:)';
sim_info.success = true;
sim_info.message = '';

rhs = @(tt, xx) local_rhs(xx, Xi, poly_order);
odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Unstable candidate models (common across a wide hyperparameter sweep)
% cause ode45 to hit finite-time blow-up, which it reports as a console
% warning rather than an error. That failure is already detected below
% (truncated/non-finite output -> sim_info.success = false), so the
% warning itself is just noise across a large sweep -- suppress it here,
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

    if isempty(X_sim) || any(~isfinite(X_sim(:))) || any(abs(X_sim(:)) > 1e6)
        sim_info.success = false;
        sim_info.message = 'Trajectory diverged (non-finite or unbounded values encountered).';
    elseif nargin >= 5 && ~isempty(t_eval) && numel(t_sim) < numel(t_eval)
        sim_info.success = false;
        sim_info.message = 'Integration terminated early (blow-up before reaching requested end time).';
    end
catch ME
    t_sim = [];
    X_sim = [];
    sim_info.success = false;
    sim_info.message = sprintf('Integration failed: %s', ME.message);
end

end


function dxdt = local_rhs(x, Xi, poly_order)
    xrow = x(:)'; % build_polynomial_library expects (N x D); N=1 for a single state
    Theta_row = build_polynomial_library(xrow, poly_order);
    dxdt = (Theta_row * Xi)';
end