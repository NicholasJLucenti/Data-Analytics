function [recommended_method, dynamics_report] = select_interpolation(t, x_channel, t_target)
    % SELECT_INTERPOLATION Adaptive interpolation choice with a
    % target-grid extrapolation guardrail.
    %
    % Inputs:
    %   t         - Original time vector (N x 1)
    %   x_channel - State variable trajectory (N x 1)
    %   t_target  - The target time vector we are mapping onto (M x 1)
    %
    % Outputs:
    %   recommended_method - one of 'linear' | 'pchip' | 'spline' | 'fourier'
    %   dynamics_report     - struct describing why that choice was made

    t = t(:); x_channel = x_channel(:);
    dynamics_report = struct();

    data_min = min(x_channel);
    data_max = max(x_channel);
    data_range = data_max - data_min;

    %% 1. Check for monotonicity / trend
    dx = diff(x_channel);
    if isempty(dx)
        dynamics_report.trend_ratio = 1;
    else
        num_positive = sum(dx > 0);
        num_negative = sum(dx < 0);
        dynamics_report.trend_ratio = max(num_positive, num_negative) / length(dx);
    end

    %% 2. Test for oscillatory dynamics (toolbox-free linear detrend)
    x_detrend = local_linear_detrend(t, x_channel);
    zero_crossings = sum(diff(x_detrend > 0) ~= 0);

    if zero_crossings >= 3 && dynamics_report.trend_ratio < 0.85
        dynamics_report.is_oscillatory = true;
    else
        dynamics_report.is_oscillatory = false;
    end

    %% 3. Initial assignment logic
    if dynamics_report.is_oscillatory
        if license('test', 'curve_fitting')
            recommended_method = 'fourier';
            dynamics_report.reason = 'Oscillatory dynamics detected.';
        else
            recommended_method = 'spline';
            dynamics_report.reason = 'Oscillatory dynamics detected (spline fallback, no Curve Fitting Toolbox).';
        end
    elseif dynamics_report.trend_ratio > 0.90
        recommended_method = 'pchip';
        dynamics_report.reason = 'Strong monotonic trend or flat baseline.';
    else
        recommended_method = 'linear';
        dynamics_report.reason = 'Standard profile. Linear selected.';
    end

    %% 4. Horizon boundary guardrail
    % Test the interpolation over the actual target timeline to catch
    % extrapolation blow-ups before they reach the rest of the pipeline.
    try
        if strcmp(recommended_method, 'fourier')
            fit_fourier = fit(t, x_channel, 'fourier1');
            test_val = fit_fourier(t_target);
        else
            test_val = interp1(t, x_channel, t_target, recommended_method, 'extrap');
        end

        if max(test_val) > (data_max + 0.5*data_range) || min(test_val) < (data_min - 0.5*data_range)
            recommended_method = 'pchip';
            dynamics_report.reason = 'CRITICAL OVERRIDE: Extrapolation exploded on the target grid. Forcing safe PCHIP.';
        end
    catch
        recommended_method = 'linear';
        dynamics_report.reason = 'Execution error during guardrail test. Falling back to linear.';
    end
end


function xd = local_linear_detrend(t, x)
    % Manual degree-1 detrend using only base MATLAB (polyfit/polyval),
    % so this file has no Signal Processing Toolbox dependency.
    t = t(:); x = x(:);
    if numel(t) < 2
        xd = x - mean(x);
        return
    end
    p = polyfit(t, x, 1);
    xd = x - polyval(p, t);
end