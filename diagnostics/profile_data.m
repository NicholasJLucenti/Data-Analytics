function report = profile_data(t, X)
    % PROFILE_DATA Master diagnostic function for SINDy Analytics.
    %
    % Inputs:
    %   t - Time vector (N x 1)
    %   X - State matrix (N x D)

    fprintf('=========================================\n');
    fprintf('        SINDY DATA ANALYTICS PROFILER    \n');
    fprintf('=========================================\n');

    [N, D] = size(X);
    report = struct();
    report.num_samples = N;
    report.dimensions = D;

    %% 1. Call Sampling Worker
    sampling = check_sampling(t);
    report.is_regularly_sampled = sampling.is_regular;
    report.dt_mean = sampling.dt_mean;
    
    fprintf('Samples: %d | State Dimensions: %d\n', N, D);
    if ~report.is_regularly_sampled
        fprintf('[WARNING] Irregular time steps detected! (CV = %.4f)\n', sampling.dt_variability);
        fprintf('          Recommendation: Run preprocessing/interpolate_missing.m\n');
    else
        fprintf('[INFO] Data is regularly sampled (dt = %.4f).\n', report.dt_mean);
    end

    %% 2. Call Noise Estimation Worker
    report.snr_db = estimate_noise(X);
    report.mean_snr = mean(report.snr_db);
    
    for i = 1:D
        fprintf('  -> State %d Estimated SNR: %.2f dB\n', i, report.snr_db(i));
    end

    %% 3. Triage & Strategy Decision Engine
    fprintf('\n--- PRELIMINARY TRIAGE REPORT ---\n');
    if report.mean_snr < 20
        report.suggested_variant = 'Weak SINDy (WSINDy)';
        fprintf('CRITICAL: Low SNR (< 20dB) detected. Bypass standard differentiation.\n');
        fprintf('STRATEGY: Use weak/integral formulation to handle noise.\n');
    elseif report.mean_snr >= 20 && report.mean_snr < 40
        report.suggested_variant = 'Standard SINDy with Regularized Differentiation';
        fprintf('STRATEGY: Use Total Variation Regularization or Savitzky-Golay for dx/dt.\n');
    else
        report.suggested_variant = 'Standard SINDy';
        fprintf('STRATEGY: Clean data detected. Finite differences + STLSQ optimizer should suffice.\n');
    end
    fprintf('=========================================\n');
end