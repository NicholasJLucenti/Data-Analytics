function period_info = estimate_trajectory_period(t, X)
%ESTIMATE_TRAJECTORY_PERIOD Estimate the dominant oscillation period of a
%trajectory, per channel and aggregated.
%
%   period_info = ESTIMATE_TRAJECTORY_PERIOD(t, X)
%
%   Reuses diagnostics/estimate_sparsity.m's FFT-peak method (the same
%   one used to measure the real data's period during diagnostics), so a
%   simulated model's period is measured identically to how the real
%   data's period was measured -- making the two directly comparable
%   rather than comparing two different period-estimation methods.
%
%   Inputs:
%     t - time vector (N x 1)
%     X - state matrix (N x D)
%
%   Output: period_info struct
%     .per_channel - 1xD vector of periods (NaN per channel if no clear
%                    dominant frequency was found -- e.g. a decaying or
%                    diverging trajectory)
%     .mean_period - mean of the finite per-channel periods (NaN if none
%                    are finite)

D = size(X, 2);
per_channel = nan(1, D);

for d = 1:D
    s = estimate_sparsity(t, X(:, d));
    per_channel(d) = s.dominant_period;
end

period_info.per_channel = per_channel;

finite_vals = per_channel(isfinite(per_channel));
if isempty(finite_vals)
    period_info.mean_period = NaN;
else
    period_info.mean_period = mean(finite_vals);
end

end