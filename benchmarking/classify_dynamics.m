function dynamics_class = classify_dynamics(t_sim, X_sim, sim_info)
%CLASSIFY_DYNAMICS Qualitatively label a simulated trajectory's long-run
%behavior, so grid-search candidates can be filtered for whether they
%preserve the qualitative dynamics (e.g. sustained oscillation) seen in
%the real data, not just whether they fit it pointwise.
%
%   dynamics_class = CLASSIFY_DYNAMICS(t_sim, X_sim, sim_info)
%
%   Inputs:
%     t_sim, X_sim - output of simulate_trajectory.m
%     sim_info     - output of simulate_trajectory.m (checked for .success)
%
%   Output: one of
%     'diverged'    - simulation failed, blew up, or is too short to judge
%     'fixed_point' - amplitude collapses toward zero in the back half of
%                     the simulation (settles to a steady state) for
%                     every channel
%     'oscillatory' - at least one channel retains substantial amplitude
%                     in the back half AND shows repeated sign changes
%                     after detrending (sustained periodic-ish motion,
%                     not just a slow monotonic drift)
%     'unknown'     - neither collapsing nor clearly oscillatory; often a
%                     slow transient that would need a longer simulation
%                     horizon to classify confidently
%
%   This is a heuristic, not a rigorous limit-cycle test (no Poincare
%   section, no Floquet analysis) -- it's meant as a cheap, automatic
%   filter across a large grid-search sweep, not a final proof of
%   dynamical structure. Treat 'oscillatory' as "worth a closer look",
%   not as a certified limit cycle.

if nargin < 3 || isempty(sim_info)
    sim_info = struct('success', true);
end

if ~sim_info.success || isempty(X_sim) || size(X_sim, 1) < 10
    dynamics_class = 'diverged';
    return
end

N = size(X_sim, 1);
splitIdx = round(N / 2);
early = X_sim(1:splitIdx, :);
late = X_sim(splitIdx+1:end, :);

earlyAmp = max(early, [], 1) - min(early, [], 1);
lateAmp = max(late, [], 1) - min(late, [], 1);

earlyAmp(earlyAmp == 0) = eps; % guard divide-by-zero for a channel with no early variation
ampRatio = lateAmp ./ earlyAmp;

if all(ampRatio < 0.1)
    dynamics_class = 'fixed_point';
    return
end

isOscillatory = false;
tLate = t_sim(splitIdx+1:end);
for d = 1:size(X_sim, 2)
    if ampRatio(d) >= 0.5
        p = polyfit(tLate, late(:, d), 1);
        detrended = late(:, d) - polyval(p, tLate);
        crossings = sum(diff(detrended > 0) ~= 0);
        if crossings >= 3
            isOscillatory = true;
            break
        end
    end
end

if isOscillatory
    dynamics_class = 'oscillatory';
else
    dynamics_class = 'unknown';
end

end