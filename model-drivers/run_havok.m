function havok = run_havok(t, y, num_delays, rank_r, varargin)
%RUN_HAVOK Hankel Alternative View of Koopman: recover a near-linear
%model in delay-coordinate space from a single (possibly partially
%observed) time series.
%
%   havok = RUN_HAVOK(t, y, num_delays, rank_r, ...)
%
%   HAVOK addresses a different failure mode than every SINDy variant
%   here: those all assume every relevant state variable is already a
%   measured channel. When that's not true -- when the real dynamics live
%   in more dimensions than you have sensors for -- SINDy has no way to
%   discover the missing dimensions; it can only fit whatever library you
%   hand it to the channels you do have. HAVOK instead builds new
%   candidate "state" coordinates directly from the time-history of a
%   single channel (via Takens delay embedding,
%   preprocessing/build_hankel_embedding.m), then finds a best-fit LINEAR
%   model in that reconstructed space:
%       d/dt v_1..r-1 = A * v_1..r-1 + B * v_r
%   where v_1..r are the leading SVD (delay-coordinate) modes and v_r,
%   the last retained mode, is treated as an external forcing term. Brunton
%   et al.'s key empirical finding is that this forcing term is often
%   near-zero most of the time with sharp intermittent bursts -- these
%   bursts are exactly where the true dynamics are nonlinear, so HAVOK
%   turns "where is this system nonlinear" into a direct, visualizable
%   question about v_r's time series, rather than requiring you to guess
%   the right nonlinear library terms up front.
%
%   Inputs:
%     t          - time vector (N x 1), uniformly sampled
%     y          - single channel (N x 1) to embed (see
%                  build_hankel_embedding.m for guidance on which channel
%                  to pick if you have more than one)
%     num_delays - number of delay copies to stack (Hankel matrix rows).
%                  A common rule of thumb is enough delays to span
%                  roughly one period of the dominant oscillation (see
%                  diagnostics/estimate_sparsity.m's dominant_period) --
%                  too few under-resolves the attractor, too many mostly
%                  adds noise and cost.
%     rank_r     - number of SVD modes to retain (including the forcing
%                  mode). rank_r-1 modes form the linear system; the
%                  r-th is the forcing term. Must be < num_delays.
%
%   Name-value options:
%     'Threshold' - fraction of the max |B| forcing magnitude below which
%                   the forcing term is treated as inactive, purely for
%                   the .forcing_active_fraction diagnostic (default 0.05)
%
%   Output: havok struct
%     .A, .B          - the fitted linear system: dV/dt(:,1:r-1) ~= V(:,1:r-1)*A' + V(:,r)*B'
%                       (A is (r-1)x(r-1), B is (r-1)x1)
%     .U, .S, .V       - the Hankel matrix's SVD (H = U*S*V'); U's columns
%                       are the delay-coordinate basis vectors, useful for
%                       reconstructing y from V
%     .V_modes         - the r retained columns of V (the mode time series
%                       used to fit A/B), each column already scaled by
%                       its singular value (i.e. V_modes = V(:,1:r)*S(1:r,1:r))
%     .t_embed         - time vector matching V_modes' rows
%     .forcing_active_fraction - fraction of samples where |v_r| exceeds
%                       Threshold * max(|v_r|) -- a quick summary of how
%                       intermittent the forcing (and therefore the
%                       "hidden nonlinearity") actually is
%     .num_delays, .rank_r - bookkeeping

p = inputParser;
addParameter(p, 'Threshold', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
parse(p, varargin{:});
opts = p.Results;

if rank_r >= num_delays
    error('run_havok:badRank', 'rank_r (%d) must be less than num_delays (%d).', rank_r, num_delays);
end
if rank_r < 2
    error('run_havok:badRank', 'rank_r must be at least 2 (1 linear mode + 1 forcing mode).');
end

%% 1. Delay embedding
[H, t_embed] = build_hankel_embedding(t, y, num_delays);

%% 2. SVD -- U's columns are the delay-coordinate basis, V's rows are the
%% mode-space trajectory
[U, S, V] = svd(H, 'econ');

r = min(rank_r, size(V, 2));
V_modes = V(:, 1:r) * S(1:r, 1:r); % scale by singular values so mode magnitudes are meaningful

%% 3. Estimate d/dt of the mode-space trajectory (finite differences on a
%% uniform grid -- t_embed inherits t's spacing from build_hankel_embedding.m)
dt = mean(diff(t_embed));
dV = zeros(size(V_modes));
for i = 1:r
    dV(:, i) = gradient(V_modes(:, i), dt);
end

%% 4. Fit the linear system dV(:,1:r-1) ~= V(:,1:r-1)*A' + V(:,r)*B'
%% i.e. [A, B] = dV(:,1:r-1)' / [V(:,1:r-1), V(:,r)]'  (least squares)
regressors = [V_modes(:, 1:r-1), V_modes(:, r)];
target = dV(:, 1:r-1);
AB = lsqminnorm(regressors, target); % (r x (r-1)), toolbox-free, robust to collinearity

A = AB(1:r-1, :)';
B = AB(r, :)';

%% 5. Forcing-intermittency diagnostic
forcing = V_modes(:, r);
maxForce = max(abs(forcing));
if maxForce > 0
    havok.forcing_active_fraction = mean(abs(forcing) > opts.Threshold * maxForce);
else
    havok.forcing_active_fraction = 0;
end

havok.A = A;
havok.B = B;
havok.U = U;
havok.S = S;
havok.V = V;
havok.V_modes = V_modes;
havok.t_embed = t_embed;
havok.num_delays = num_delays;
havok.rank_r = r;

end