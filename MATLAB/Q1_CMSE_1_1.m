%% Q1 Classical and Modern Spectrum Estimation
%{
Show analytically and through simulations that the definition of PSD in (7) 
is equivalent to that in (9) under a mild [5] assumption that the 
covariance sequence r(k) decays rapidly, that is [1]

Provide a simulation for the case when this equivalence does not hold.
 Explain the reasons.
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 1;
q_initialise;

SAVE_FIGS = true;

%% 
% Set rng seed for consistency
rSeed = 13;
rng(rSeed)

% Provide a simulation for the case when this equivalence does not hold.
% Explain the reasons.
fs = 50; % Sampling frequency of both simulations
N = 500; % Select initial value of N
t = 0:1/fs:(N-1)/fs;
f_wave = 2; % Frequency of wave
wave = sin(2*pi*f_wave*t) + 2*randn(1,N);
% wave

% Find autocorrelation
R_wave_hat  = xcorr(wave,'unbiased'); % Find unbiased estimate
R_wave_hatB = xcorr(wave,'biased');   % Find biased estimate
maxlag = length(t) - 1; % Find Maximum lag
tau = -maxlag:maxlag;

% find PSD
fft_power_sf  = fft(wave).*conj(fft(wave))/N;    % signal
fft_power_cr  = fft(R_wave_hat);                 % unbiased
fft_power_crB = fft(R_wave_hatB);                % biased

f = 0:fs/N:fs-1/N;

%% Plots
close all;
fH = []; % reset the figure handle array
figure;
    stem(fft_power_cr);

fH{length(fH) +1} = figure; hold on
    stem(tau, R_wave_hat, 'DisplayName', 'Unbiased');
    stem(tau, R_wave_hatB, 'DisplayName', 'Biased');
    % axis([-40 40 -inf inf]);

    title(sprintf('ACF of a Sine Wave \n $ sin(2\\pi(%i)t) + 2(\\eta(t))$',f_wave));
    xlabel('Lags');
    ylabel('Magnitude', 'Color','k');
    grid minor;
    legend('show', 'Location','best')
    
% Plot PSD for N=500    
fH{length(fH) +1} = figure; hold on
    stem(f,fft_power_cr(N:end), 'DisplayName', 'Unbiased');
    stem(f,fft_power_crB(N:end), 'DisplayName','Biased');
    stem(f,fft_power_sf, 'DisplayName','Actual');
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    title(sprintf('PSD Comparison for ACF Estimators\n $ sin(2\\pi(%i)t) + 2(\\eta(t))$',f_wave))
    grid minor;
    legend('show', 'Location','best')
    

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_1b_fig',num2str(ii,'%02i')],'pdf')
    end
end