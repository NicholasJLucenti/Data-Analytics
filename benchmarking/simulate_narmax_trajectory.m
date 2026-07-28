function [t_sim, X_sim, sim_info] = simulate_narmax_trajectory(Xi, lag_order, library_spec, t, X_seed, varargin)
%SIMULATE_NARMAX_TRAJECTORY Forward-roll a fitted NARMAX model step by
%step across a time grid.
%
%   [t_sim, X_sim, sim_info] = SIMULATE_NARMAX_TRAJECTORY(Xi, lag_order, library_spec, t, X_seed, ...)
%
%   Unlike variants/run_standard_sindy.m or run_weak_sindy.m's models,
%   NARMAX has no notion of a single initial condition -- it needs
%   lag_order PAST states to predict the next one. This rolls the fitted
%   discrete-time map forward: X_k = Theta([X_{k-1},...,X_{k-lag_order}]) * Xi,
%   using its OWN previous predictions as the lag history (not the real
%   data), so this is a genuine forward simulation and not one-step-ahead
%   prediction.
%
%   IMPORTANT: the lag-history column ordering here (lag varies slower
%   than channel: column = (lag-1)*D + d) must exactly match
%   variants/run_narmax.m's Z construction, or Theta*Xi will be applied
%   to a differently-ordered regressor than the one Xi was fit against.
%
%   Inputs:
%     Xi           - coefficient matrix from run_narmax.m (M x D)
%     lag_order    - the lag_order run_narmax.m was fit with
%     library_spec - the exact library_spec run_narmax.m was fit with
%     t            - full target time vector (N x 1) to roll forward across
%     X_seed       - (lag_order x D) matrix of the first lag_order REAL
%                    states, used to seed the rollout
%
%   Name-value options:
%     'DivergenceBound' - stop early if any state exceeds this magnitude
%                         (default 1e6, matches simulate_trajectory.m)
%     'MaxWallSeconds'  - hard wall-clock cap on the whole rollout
%                         (default 3, matches simulate_trajectory.m)
%
%   Outputs: t_sim, X_sim, sim_info -- same contract as
%   benchmarking/simulate_trajectory.m (sim_info.success/.message; X_sim
%   is truncated and success=false if the rollout diverges or times out),
%   so compute_narmax_trajectory_error.m, classify_dynamics.m, and
%   estimate_trajectory_period.m all work on NARMAX output unchanged.

p = inputParser;
addParameter(p, 'DivergenceBound', 1e6, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MaxWallSeconds', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});
opts = p.Results;

t = t(:);
N = numel(t);
D = size(X_seed, 2);

if size(X_seed, 1) ~= lag_order
    error('simulate_narmax_trajectory:badSeed', 'X_seed must have exactly lag_order (%d) rows.', lag_order);
end
if N <= lag_order
    error('simulate_narmax_trajectory:tooFewPoints', 't must have more than lag_order points.');
end

X_sim = zeros(N, D);
X_sim(1:lag_order, :) = X_seed;

sim_info.success = true;
sim_info.message = '';

startTime = tic;

for k = (lag_order + 1):N
    if toc(startTime) > opts.MaxWallSeconds
        sim_info.success = false;
        sim_info.message = sprintf('NARMAX rollout exceeded wall-clock budget (%.1fs).', opts.MaxWallSeconds);
        t_sim = t(1:k-1);
        X_sim = X_sim(1:k-1, :);
        return
    end

    zrow = zeros(1, D * lag_order);
    col = 0;
    for lag = 1:lag_order
        for d = 1:D
            col = col + 1;
            zrow(col) = X_sim(k - lag, d);
        end
    end

    Theta_row = build_library(zrow, library_spec);
    X_sim(k, :) = Theta_row * Xi;

    if any(~isfinite(X_sim(k, :))) || any(abs(X_sim(k, :)) > opts.DivergenceBound)
        sim_info.success = false;
        sim_info.message = 'NARMAX rollout diverged (non-finite or unbounded state).';
        t_sim = t(1:k);
        X_sim = X_sim(1:k, :);
        return
    end
end

t_sim = t;

end