function sampling_report = check_sampling(t)
    % CHECK_SAMPLING Analyzes the time vector for step uniformity.
    %
    % Input:
    %   t - Vector of time stamps (N x 1)
    %
    % Output:
    %   sampling_report - Struct containing timing metrics

    sampling_report = struct();
    
    dt = diff(t);
    sampling_report.dt_mean = mean(dt);
    sampling_report.dt_std = std(dt);
    
    % Coefficient of Variation (CV) for time steps
    sampling_report.dt_variability = sampling_report.dt_std / sampling_report.dt_mean;
    
    % Strict threshold for uniformity (allowing minor floating-point errors)
    if sampling_report.dt_variability > 1e-4
        sampling_report.is_regular = false;
    else
        sampling_report.is_regular = true;
    end
end