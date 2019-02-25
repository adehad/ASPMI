%% Q1 Classical and Modern Spectrum Estimation
%  1.2a
%{
    Apply one periodogram-based spectral estimation technique
(possibly after some preprocessing) to the sunspot time [10] series. 
    Explain what aspect of the spectral estimate changes when the
mean and trend from the data are removed (use the MATLAB commands mean and detrend).
    Explain how the perception of the periodicities in the data changes when the
data is transformed by:
 ?rst applying the logarithm to each data sample and
 then subtracting the sample mean from this logarithmic data.
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;
% 
% % Initialise common functions for this question
questionNum = 1;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
% SAVE_FIGS = true;


%% Initialising Input Data and Variables

sunspot = load('sunspot.dat');
% dim 1: years (increments in 1 year elements)
% dim 2: Zurich Relative Sunspot Number (related to number and size of spots)

years = sunspot(:,1);
relZNums = sunspot(:,2);

% Recenter lambda function (removes mean bias)
recenter = @(x) (x - mean(x));

% Preprocess data
relZNums_prePro = detrend(recenter(relZNums));

% Preprocess data (natural log)
relZNums_log = log(relZNums +eps); % +eps prevent log(0)
relZNums_log_recenter = recenter(relZNums);

fs = 1; % 1 year per sample
N = size(sunspot,1); % Total Signal Length
K = 2^12; % 4096 - a zero-padded signal size
    % Why we zero pad can also be seen here:
    % https://dsp.stackexchange.com/questions/741/why-should-i-zero-pad-a-signal-before-taking-the-fourier-transform
%% Initialising a Window
% Requires Signal Processing ToolBox
% https://uk.mathworks.com/help/signal/windows.html

% Run: windowDesigner
% for an interactive look at window types

win_BKH = blackmanharris(N); % Blackman Harris Window
winLabel = ["BK-Harris"];

% win_

%% Periodogram / PSD Calculations
legendString = [];

% Raw Sunspot Data
FFT_sunspot = (relZNums .* win_BKH)'; % apply the window
FFT_sunspot = [ FFT_sunspot zeros(1, K-N) ]; % zero padd the signal
FFT_sunspot = fftshift(fft(FFT_sunspot)); % apply the fft
P_sunspot = pow2db( abs(FFT_sunspot).^2 /(N *2*pi) )'; % convert to power spectrum (db)
legendString = [legendString, "Raw SunSpots"];

% Detrended+Recentered Sunspot Data
FFT_sunspot_prePro = (relZNums_prePro .* win_BKH)'; % apply the window
FFT_sunspot_prePro = [ FFT_sunspot_prePro zeros(1, K-N) ]; % zero padd the signal
FFT_sunspot_prePro = fftshift(fft(FFT_sunspot_prePro)); % apply the fft
P_sunspot_prePro = pow2db( abs(FFT_sunspot_prePro).^2 /(N *2*pi) )'; % convert to power spectrum (db)
legendString = [legendString, "$-\mathtt{mean}$ \& $\mathtt{detrend}$"];

% Log Sunspot Data
FFT_sunspot_log = (relZNums_log .* win_BKH)'; % apply the window
FFT_sunspot_log = [ FFT_sunspot_log zeros(1, K-N) ]; % zero padd the signal
FFT_sunspot_log = fftshift(fft(FFT_sunspot_log)); % apply the fft
P_sunspot_log = pow2db( abs(FFT_sunspot_log).^2 /(N *2*pi) )'; % convert to power spectrum (db)
legendString = [legendString, "$\mathtt{log}$"];

% Log+Recentered Sunspot Data
FFT_sunspot_log_recenter = (relZNums_log_recenter .* win_BKH)'; % apply the window
FFT_sunspot_log_recenter = [ FFT_sunspot_log_recenter zeros(1, K-N) ]; % zero padd the signal
FFT_sunspot_log_recenter = fftshift(fft(FFT_sunspot_log_recenter)); % apply the fft
P_sunspot_log_recenter = pow2db( abs(FFT_sunspot_log_recenter).^2 /(N *2*pi) )'; % convert to power spectrum (db)
legendString = [legendString, "$\mathtt{log}$ \& $-\mathtt{mean}$"];


% Establish a logical array to generate one sided PSD
oneSide = [zeros(1,K/2) ones(1,K/2)]; % all spectrums the same length
oneSide = (oneSide == 1); % converts to boolean

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable

fAx = 0:2*fs/K:1-fs/K; % normalised frequency axis

fH{1} = figure;
    hold on
    plot(fAx, P_sunspot(oneSide)); 
%     drawnow; pause;
    plot(fAx, P_sunspot_prePro(oneSide));
%     drawnow; pause;
%     plot(fAx, P_sunspot_log(oneSide));
%     drawnow; pause;
    plot(fAx, P_sunspot_log_recenter(oneSide));
%     drawnow; pause;
    title(sprintf("{%s} Window Periodogram ", winLabel));
    xlabel("Normalised Frequency ($\frac{\pi\ rad}{sample}$)");
    ylabel("Power Density (dB)");
    legend(legendString([1,2,4]))
    grid minor
    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_2a_fig',num2str(ii,'%02i')],'pdf')
    end
end

