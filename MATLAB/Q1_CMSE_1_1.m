%% Q1 Classical and Modern Spectrum Estimation
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 1;
q_initialise;

SAVE_FIGS = true;

%% 
% Provide a simulation for the case when this equivalence does not hold.
% Explain the reasons.
fs = 50; % Sampling frequency of both simulations
N = 500; % Select initial value of N
t = 0:1/fs:(N-1)/fs;
f_wave = 2; % Frequency of wave
wave = sin(2*pi*f_wave*t) + 2*randn(1,N);
% wave

% Find autocorrelation
R_wave_hat = xcorr(wave,'unbiased'); % Find unbiased estimate
maxlag = length(t) - 1; % Find Maximum lag
tau = -maxlag:maxlag;
figure;
stem(tau, R_wave_hat);
% axis([-40 40 -inf inf]);
title('ACF of a Sine Wave');
xlabel('Lags');
ylabel('Autocorrelation');

% Plot PSD for N=500
N = 500;
t = 0:1/fs:(N-1)/fs;
% wave = sin(2*pi*f_wave*t) + 2*randn(1,N);
fft_power_cr = fft(R_wave_hat);
fft_power_sf = fft(wave).*conj(fft(wave))/N;
f = 0:fs/N:fs-1/N;
fH{1} = figure; hold on
stem(f,fft_power_sf);
stem(f,fft_power_cr(end/2:end));
xlabel('Frequency (Hz)')
ylabel('PSD')
title('Limiting Case: Unbiased ACF Estimators')

figure;
stem(fft_power_cr);

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_1b_fig',num2str(ii,'%02i')],'pdf')
    end
end