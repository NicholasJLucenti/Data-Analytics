function snr_db = estimate_noise(X)
    % ESTIMATE_NOISE Estimates Signal-to-Noise Ratio for multi-dimensional data.
    %
    % Inputs:
    %   X - Matrix of state trajectories (N x D)
    %
    % Output:
    %   snr_db - Vector of SNR values in dB for each dimension (D x 1)

    [N, D] = size(X);
    snr_db = zeros(D, 1);
    
    % Guardrail: If dataset is too short, FFT-based noise profiling fails
    if N < 16
        warning('Dataset is too short (N = %d) for spectral noise estimation. Defaulting to variance-based heuristic.', N);
        for i = 1:D
            signal = X(:, i);
            % Use a basic high-pass finite difference acting as a rough noise proxy
            noise_proxy = diff(signal);
            signal_power = var(signal);
            noise_power = var(noise_proxy);
            if noise_power > 0
                snr_db(i) = 10 * log10(signal_power / noise_power);
            else
                snr_db(i) = Inf;
            end
        end
        return;
    end

    %% Spectral Estimation for Sufficiently Long Data
    for i = 1:D
        signal_detrend = detrend(X(:, i));
        fft_res = fft(signal_detrend);
        psd = abs(fft_res(1:floor(N/2)+1)).^2;
        
        len_psd = length(psd);
        
        % Dynamic indexing to ensure we never grab empty arrays
        signal_idx = 1:max(1, floor(len_psd * 0.15));
        noise_idx = min(len_psd, ceil(len_psd * 0.85)):len_psd;
        
        signal_power = mean(psd(signal_idx));
        noise_power = mean(psd(noise_idx));
        
        if noise_power > 0 && ~isnan(signal_power)
            snr_db(i) = 10 * log10(signal_power / noise_power);
        else
            snr_db(i) = Inf; 
        end
    end
end