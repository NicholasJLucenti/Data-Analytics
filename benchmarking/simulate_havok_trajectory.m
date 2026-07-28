function [t_sim, y_sim, sim_info] = simulate_havok_trajectory(havok, t_eval)
%SIMULATE_HAVOK_TRAJECTORY Forward-integrate a fitted HAVOK model and
%reconstruct the original (single-channel) trajectory from it.
%
%   [t_sim, y_sim, sim_info] = SIMULATE_HAVOK_TRAJECTORY(havok, t_eval)
%
%   Integrates dV/dt(1:r-1) = A*V(1:r-1) + B*v_r forward, using the
%   model's OWN fitted forcing history v_r(t) (interpolated onto whatever
%   times the integrator asks for) rather than re-deriving it -- HAVOK's
%   forcing term isn't itself modeled by the linear system, it's the
%   piece treated as external input, so reusing its recorded time series
%   is the standard way to validate the fit (this is a reconstruction
%   check, not a fully autonomous forecast beyond the fitted forcing
%   record).
%
%   Inputs:
%     havok  - struct from run_havok.m
%     t_eval - times to evaluate at (must lie within havok.t_embed's
%              range, since the forcing term is only defined there)
%
%   Outputs:
%     t_sim, sim_info - same success/message contract as
%                        benchmarking/simulate_trajectory.m
%     y_sim            - reconstructed single-channel trajectory (numel(t_eval) x 1),
%                        obtained by integrating the mode-space system
%                        then projecting the leading mode back through
%                        havok.U's first row scaling (the standard HAVOK
%                        reconstruction: y ~= U(1,1:r-1) * V(:,1:r-1)')

t_eval = t_eval(:);
r = havok.rank_r;

if any(t_eval < havok.t_embed(1) - eps) || any(t_eval > havok.t_embed(end) + eps)
    warning('simulate_havok_trajectory:extrapolating', ...
        'Some requested times fall outside the forcing term''s fitted range; forcing will be held at the nearest edge value.');
end

forcing_interp = @(tt) interp1(havok.t_embed, havok.V_modes(:, r), tt, 'linear', 'extrap');

v0 = interp1(havok.t_embed, havok.V_modes(:, 1:r-1), t_eval(1), 'linear', 'extrap')';

rhs = @(tt, vv) havok.A * vv + havok.B * forcing_interp(tt);

sim_info.success = true;
sim_info.message = '';

try
    odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_sim, V_sim] = ode45(rhs, t_eval, v0, odeOpts);

    if isempty(V_sim) || any(~isfinite(V_sim(:))) || size(V_sim, 1) ~= numel(t_eval)
        sim_info.success = false;
        sim_info.message = 'HAVOK reconstruction diverged or terminated early.';
        y_sim = [];
        return
    end

    % Standard HAVOK reconstruction: project the leading delay-coordinate
    % modes back through U's first row (the "current time" row of the
    % delay embedding) to recover the original channel's scale.
    y_sim = V_sim * havok.U(1, 1:r-1)';

catch ME
    t_sim = [];
    y_sim = [];
    sim_info.success = false;
    sim_info.message = sprintf('HAVOK integration failed: %s', ME.message);
end

end