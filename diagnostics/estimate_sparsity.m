function sparsity = estimate_sparsity(t, y)
%ESTIMATE_SPARSITY Estimate how sparsely a signal is sampled relative to
%its own characteristic dynamical timescale.
%
%   sparsity = ESTIMATE_SPARSITY(t, y) returns a struct:
%     .n_points           - number of unique finite samples
%     .dominant_period    - estimated period of the dominant oscillatory/
%                           trend component (via FFT peak), in the same
%                           units as t. NaN if no clear dominant frequency
%                           could be identified.
%     .points_per_period  - average number of samples covering one cycle
%                           of the dominant dynamics
%     .is_sparse          - true if points_per_period falls below a
%                           minimum-safe threshold (default 10), or if
%                           there are too few points to tell
%
%   Rationale: 16 points over a 12-hour window is sparse if the dynamics
%   oscillate every 2 hours (6 points/cycle) but dense if the dynamics
%   are monotonic and slow over the same window. This metric is what
%   actually determines whether derivative estimation / densification is
%   safe, rather than raw point count alone.

t = t(:); y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid); y = y(valid);
[t, sortIdx] = sort(t); y = y(sortIdx);

[tu, ~, ic] = unique(t);
if numel(tu) < numel(t)
    y = accumarray(ic, y, [], @mean);
    t = tu;
end

sparsity.n_points = numel(t);

if numel(t) < 8 || t(end) == t(1)
    sparsity.dominant_period = NaN;
    sparsity.points_per_period = NaN;
    sparsity.is_sparse = true; % too few points / no span to judge otherwise -> assume worst case
    return
end

nGrid = max(64, 2^nextpow2(numel(t)));
tg = linspace(t(1), t(end), nGrid)';
yg = interp1(t, y, tg, 'linear');
yg = yg - mean(yg);

if var(yg) == 0
    sparsity.dominant_period = NaN;
    sparsity.points_per_period = NaN;
    sparsity.is_sparse = false; % constant signal has no meaningful timescale to be sparse relative to
    return
end

Y = fft(yg);
P = abs(Y).^2;
halfN = floor(nGrid/2);
P = P(2:halfN); % drop DC term

fs = 1 / (tg(2) - tg(1));
freqs = fs * (1:(halfN-1))' / nGrid;

[~, peakIdx] = max(P);
peakFreq = freqs(peakIdx);

if peakFreq <= 0
    sparsity.dominant_period = NaN;
    sparsity.points_per_period = NaN;
    sparsity.is_sparse = sparsity.n_points < 20; % fallback heuristic
    return
end

sparsity.dominant_period = 1 / peakFreq;

meanRawDt = (t(end) - t(1)) / (numel(t) - 1);
sparsity.points_per_period = sparsity.dominant_period / meanRawDt;

minSafePointsPerPeriod = 10;
sparsity.is_sparse = sparsity.points_per_period < minSafePointsPerPeriod;

end