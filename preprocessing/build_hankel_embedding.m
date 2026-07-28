function [H, t_embed] = build_hankel_embedding(t, y, num_delays)
%BUILD_HANKEL_EMBEDDING Construct a Hankel (delay-coordinate) matrix from
%a single time series, per Takens' embedding theorem.
%
%   [H, t_embed] = BUILD_HANKEL_EMBEDDING(t, y, num_delays)
%
%   Takens' theorem: for a sufficiently rich single time series, stacking
%   enough time-delayed copies of it reconstructs a state space
%   diffeomorphic to the true (possibly higher-dimensional, possibly
%   partially unobserved) underlying dynamics. This is the standard fix
%   for "not all relevant state variables are included in the data" --
%   rather than requiring new measurements, it recovers missing
%   dimensionality from the time-history of what you DO have.
%
%   Inputs:
%     t          - time vector (N x 1), assumed uniformly sampled (use
%                  preprocessing/align_and_truncate.m's output)
%     y          - single channel (N x 1) to embed. HAVOK is normally
%                  applied to ONE representative channel -- pick
%                  whichever is best-sampled/least noisy if you have a
%                  choice, since delay coordinates of a noisy channel
%                  compound that noise across every delay.
%     num_delays - number of delay copies to stack (the Hankel matrix has
%                  this many rows). Larger num_delays captures longer
%                  memory / can resolve more hidden dimensions, but
%                  shrinks the usable trajectory length (N - num_delays + 1)
%                  and amplifies noise sensitivity.
%
%   Outputs:
%     H       - Hankel matrix (num_delays x (N - num_delays + 1)). Row i
%               is y shifted by (i-1) samples: H(i,:) = y(i : i+cols-1)'.
%               Column j is one delay-embedded "state" at the time
%               t_embed(j).
%     t_embed - time vector (1 x (N - num_delays + 1)) matching H's
%               columns -- these are the times of the FIRST (least
%               delayed) row's samples, i.e. the "current" time for each
%               embedded state.

t = t(:);
y = y(:);

if numel(t) ~= numel(y)
    error('build_hankel_embedding:sizeMismatch', 'numel(t) must equal numel(y).');
end

N = numel(y);
if num_delays < 2
    error('build_hankel_embedding:badNumDelays', 'num_delays must be at least 2.');
end
if num_delays >= N
    error('build_hankel_embedding:tooManyDelays', ...
        'num_delays (%d) must be less than the number of samples (%d).', num_delays, N);
end

numCols = N - num_delays + 1;
H = zeros(num_delays, numCols);
for i = 1:num_delays
    H(i, :) = y(i:(i + numCols - 1))';
end

t_embed = t(1:numCols)';

end