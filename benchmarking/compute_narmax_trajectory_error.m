function metrics = compute_narmax_trajectory_error(t, X, Xi, lag_order, library_spec)
%COMPUTE_NARMAX_TRAJECTORY_ERROR Forward-roll a fitted NARMAX model from
%the data's own seed history and score it against the real trajectory.
%
%   metrics = COMPUTE_NARMAX_TRAJECTORY_ERROR(t, X, Xi, lag_order, library_spec)
%
%   Mirrors benchmarking/compute_trajectory_error.m's contract exactly
%   (same output fields), but calls simulate_narmax_trajectory.m's
%   discrete-time rollout instead of simulate_trajectory.m's ode45
%   integration -- NARMAX has no continuous-time ODE to integrate.
%
%   Output: metrics struct with .success, .message, .rmse,
%   .rmse_per_channel, .normalized_rmse -- identical shape to
%   compute_trajectory_error.m's output.

X_seed = X(1:lag_order, :);
[t_sim, X_sim, sim_info] = simulate_narmax_trajectory(Xi, lag_order, library_spec, t, X_seed);

metrics.success = sim_info.success;
metrics.message = sim_info.message;

if ~sim_info.success || isempty(X_sim) || size(X_sim, 1) ~= numel(t)
    metrics.rmse = Inf;
    metrics.rmse_per_channel = Inf(1, size(X, 2));
    metrics.normalized_rmse = Inf;
    return
end

diffX = X_sim - X;
metrics.rmse_per_channel = sqrt(mean(diffX.^2, 1));
metrics.rmse = sqrt(mean(diffX(:).^2));

channelStd = std(X, 0, 1);
channelStd(channelStd == 0) = 1;
metrics.normalized_rmse = mean(metrics.rmse_per_channel ./ channelStd);

end