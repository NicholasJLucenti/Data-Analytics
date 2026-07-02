function [recommended_method, dynamics_report] = select_interpolation(t, x_channel, t_target)
    % SELECT_INTERPOLATION Adaptive selection with multi-grid horizon validation.
    %
    % Inputs:
    %   t         - Original time vector (N x 1)
    %   x_channel - State variable trajectory (N x 1)
    %   t_target  - The target time vector we are mapping onto (M x 1)

    N = length(x_channel);
    dynamics_report = struct();
    
    data_min = min(x_channel);
    data_max = max(x_channel);
    data_range = data_max - data_min;

    %% 1. Check for Monotonicity / Trend 
    dx = diff(x_channel);
    num_positive = sum(dx > 0);
    num_negative = sum(dx < 0);
    dynamics_report.trend_ratio = max(num_positive, num_negative) / length(dx);
    
    %% 2. Test for Oscillatory Dynamics
    x_detrend = detrend(x_channel);
    zero_crossings = sum(diff(x_detrend > 0) ~= 0);
    
    if zero_crossings >= 3 && dynamics_report.trend_ratio < 0.85
        dynamics_report.is_oscillatory = true;
    else
        dynamics_report.is_oscillatory = false;
    end

    %% 3. Initial Assignment Logic
    if dynamics_report.is_oscillatory
        if license('test', 'curve_fitting')
            recommended_method = 'fourier'; 
            dynamics_report.reason = 'Oscillatory dynamics detected.';
        else
            recommended_method = 'spline'; 
            dynamics_report.reason = 'Oscillatory dynamics detected (spline fallback).';
        end
    elseif dynamics_report.trend_ratio > 0.90
        recommended_method = 'pchip';
        dynamics_report.reason = 'Strong monotonic trend or flat baseline.';
    else
        recommended_method = 'linear';
        dynamics_report.reason = 'Standard profile. Linear selected.';
    end

    %% 4. FIXED: Horizon Boundary Guardrail
    % Test the interpolation over the actual *target* timeline to catch extrapolation errors
    try
        if strcmp(recommended_method, 'fourier')
            fit_fourier = fit(t, x_channel, 'fourier1');
            test_val = fit_fourier(t_target);
        else
            test_val = interp1(t, x_channel, t_target, recommended_method, 'extrap');
        end
        
        % If it shoots past the raw limits anywhere on the target timeline, trigger override
        if max(test_val) > (data_max + 0.5*data_range) || min(test_val) < (data_min - 0.5*data_range)
            recommended_method = 'pchip'; 
            dynamics_report.reason = 'CRITICAL OVERRIDE: Extrapolation exploded on the target grid. Forcing safe PCHIP.';
        end
    catch
        recommended_method = 'linear';
        dynamics_report.reason = 'Execution error. Falling back to linear.';
    end
end