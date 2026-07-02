function [recommended_method, dynamics_report] = select_interpolation(t, x_channel)
    % SELECT_INTERPOLATION Analyzes a single state trajectory channel to 
    % automatically choose the safest and most accurate interpolation method.
    %
    % Inputs:
    %   t         - Time vector (N x 1)
    %   x_channel - Single column vector of a state variable trajectory (N x 1)
    %
    % Outputs:
    %   recommended_method - String ('pchip', 'spline', 'linear', or 'harmonic')
    %   dynamics_report    - Struct containing the underlying analytics

    N = length(x_channel);
    dynamics_report = struct();
    
    %% 1. Check for Monotonicity / Trend 
    % If data is strictly increasing or decreasing, polynomials will overshoot.
    dx = diff(x_channel);
    num_positive = sum(dx > 0);
    num_negative = sum(dx < 0);
    
    % If the data moves in one direction more than 85% of the time, it has a strong trend
    dynamics_report.trend_ratio = max(num_positive, num_negative) / length(dx);
    
    %% 2. Test for Oscillatory Dynamics via Zero-Crossings of Detrended Signal
    x_detrend = detrend(x_channel);
    
    % Count how many times the signal crosses its mean value
    zero_crossings = sum(diff(x_detrend > 0) ~= 0);
    dynamics_report.zero_crossings = zero_crossings;
    
    % Heuristic: If it crosses its mean frequently relative to dataset length, 
    % and doesn't have a strict monotonic trend, it's likely shifting/oscillating.
    if zero_crossings >= 3 && dynamics_report.trend_ratio < 0.85
        dynamics_report.is_oscillatory = true;
    else
        dynamics_report.is_oscillatory = false;
    end

    %% 3. Automated Selection Logic
    if dynamics_report.is_oscillatory
        % For oscillating data (like your Hes1 mRNA), splines capture peaks well.
        % If data is highly regular, a harmonic/sinusoidal approach works, but
        % 'spline' is the safest mathematical default for generalized oscillations.
        recommended_method = 'spline'; 
        dynamics_report.reason = 'Oscillatory behavior detected. Splines preserve smooth peak curvatures.';
        
    elseif dynamics_report.trend_ratio > 0.90
        % If it's a flat line or sudden explosion (like the flat protein profile),
        % use PCHIP to prevent Runge's phenomenon overshoot.
        recommended_method = 'pchip';
        dynamics_report.reason = 'Strong monotonic trend or flat baseline detected. PCHIP prevents overshoot.';
        
    else
        % Default fallback for sparse, highly irregular or uncertain profiles
        recommended_method = 'linear';
        dynamics_report.reason = 'Standard dynamic profile. Linear interpolation chosen to minimize structural assumptions.';
    end
end