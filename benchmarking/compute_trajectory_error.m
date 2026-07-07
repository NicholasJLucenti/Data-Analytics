function metrics = compute_trajectory_error(t, X, Xi, poly_order)
%COMPUTE_TRAJECTORY_ERROR Forward-simulate a discovered model from the
%data's own initial condition and score it against the real trajectory.
%
%   metrics = COMPUTE_TRAJECTORY_ERROR(t, X, Xi, poly_order)
%
%   Inputs:
%     t          - real time vector (N x 1), dense/uniform (from
%                  preprocessing/align_and_truncate.m)
%     X          - real state matrix (N x D)
%     Xi         - discovered coefficient matrix (M x D)
%     poly_order - polynomial order the library Xi was fit against
%
%   Output: metrics struct
%     .success           - whether the forward simulation stayed finite/bounded
%     .message           - explanation if not
%     .rmse              - root-mean-square error, pooled across all channels
%     .rmse_per_channel  - 1 x D
%     .normalized_rmse   - rmse_per_channel divided by that channel's own
%                          std(X), then averaged -- a scale-free score so
%                          channels with very different magnitudes (e.g.
%                          mRNA vs protein counts) contribute comparably
%
%   This is the check that in-sample fit metrics (AIC/BIC, ensemble
%   inclusion probability) cannot provide: it asks whether the discovered
%   equations actually reproduce the observed dynamics when integrated
%   forward, not just whether they locally match derivative estimates.

x0 = X(1, :);
[t_sim, X_sim, sim_info] = simulate_trajectory(Xi, poly_order, [t(1), t(end)], x0, t);

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
channelStd(channelStd == 0) = 1; % guard divide-by-zero for a constant channel
metrics.normalized_rmse = mean(metrics.rmse_per_channel ./ channelStd);

end