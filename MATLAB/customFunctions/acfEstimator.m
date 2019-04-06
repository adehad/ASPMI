function [r_biased, r_unbiased, lags, PSD_biased, PSD_unbiased, fs] = acfEstimator(x)
% Auto Correlation Estimator
% Input: 
%       - x
% Output: 
%       * r_biased = Biased Estimate
%       * r_unbiased = Unbiased_estimate
%       * lag = ,
%       * PSD_biased = Biased Periodogram/PSD, PSD_Unbiased = Unbiased
%       * fs = Frequency

    % Auto-Correlation Estimates
    [r_biased, lags] = xcorr(x, 'biased');
    [r_unbiased, ~]  = xcorr(x, 'unbiased');
    % Shift/Center ACF
    acf_biased   = ifftshift(r_biased);
    acf_unbiased = ifftshift(r_unbiased);
    % PSD using FFT
    PSD_biased   = real(fftshift(fft(acf_biased)))   /(2*pi);
    PSD_unbiased = real(fftshift(fft(acf_unbiased))) /(2*pi);
    % Normalise frequency
    fs = lags/max(lags);
    
end