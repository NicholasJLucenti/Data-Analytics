function noise = estimate_noise(t, y)
%ESTIMATE_NOISE Estimate the noise floor and SNR of a single (possibly
%unevenly sampled) channel.
%
%   noise = ESTIMATE_NOISE(t, y) returns a struct:
%     .snr_db        - estimated signal-to-noise ratio, in dB
%     .noise_power   - estimated noise variance
%     .signal_power  - estimated signal variance
%     .method        - string describing which estimator was used
%
%   Two regimes:
%     - Fewer than 16 finite points: spectral estimation is unreliable,
%       so a variance-based heuristic is used instead (noise proxy =
%       first-difference variance).
%     - 16+ points: signal is resampled onto a uniform grid, linearly
%       detrended (no toolbox required), and its power spectrum is split
%       into a low-frequency band (assumed dynamics-dominated) and a
%       high-frequency band (assumed noise-dominated).
%
%   This is a lightweight heuristic for triage/routing decisions (see
%   selection/select_route.m), not a substitute for a proper noise model.

t = t(:); y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid); y = y(valid);

[t, sortIdx] = sort(t);
y = y(sortIdx);

% collapse exact duplicate timestamps by averaging
[tu, ~, ic] = unique(t);
if numel(tu) < numel(t)
    y = accumarray(ic, y, [], @mean);
    t = tu;
end

n = numel(t);

if n < 16
    % Variance-based fallback: treat the first-difference series as a
    % rough noise proxy (differencing suppresses slow trend, leaving
    % high-frequency jitter).
    signal_power = var(y);
    if n >= 2
        noise_proxy = diff(y);
        noise_power = var(noise_proxy);
    else
        noise_power = 0;
    end

    if noise_power > 0
        noise.snr_db = 10*log10(signal_power / noise_power);
    else
        noise.snr_db = Inf;
    end
    noise.noise_power = noise_power;
    noise.signal_power = signal_power;
    noise.method = sprintf('variance-based fallback (n=%d < 16, spectral estimate unreliable)', n);
    return
end

if t(end) == t(1)
    noise.snr_db = NaN;
    noise.noise_power = NaN;
    noise.signal_power = var(y);
    noise.method = 'zero time span, cannot resample for spectral estimate';
    return
end

% resample onto a uniform grid for FFT
nGrid = max(64, 2^nextpow2(n));
tg = linspace(t(1), t(end), nGrid)';
yg = interp1(t, y, tg, 'linear');

% manual linear detrend (avoids Signal Processing Toolbox's detrend())
p = polyfit(tg, yg, 1);
yg = yg - polyval(p, tg);

total_power = var(yg);
if total_power == 0
    noise.snr_db = NaN;
    noise.noise_power = 0;
    noise.signal_power = 0;
    noise.method = 'signal has zero variance';
    return
end

Y = fft(yg);
halfN = floor(nGrid/2) + 1;
psd = abs(Y(1:halfN)).^2;
len_psd = length(psd);

% low 15% of the spectrum = dynamics band, top 15% = noise band
signal_idx = 1:max(1, floor(len_psd * 0.15));
noise_idx = min(len_psd, ceil(len_psd * 0.85)):len_psd;

signal_power = mean(psd(signal_idx));
noise_power = mean(psd(noise_idx));

if noise_power > 0 && isfinite(signal_power)
    noise.snr_db = 10*log10(signal_power / noise_power);
else
    noise.snr_db = Inf;
end
noise.noise_power = noise_power;
noise.signal_power = signal_power;
noise.method = 'low/high-band spectral power ratio (linear-detrended, uniform-resampled)';

end