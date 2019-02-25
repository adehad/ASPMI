%% Q1 Classical and Modern Spectrum Estimation
% 1.3a
%{
 Write a MATLAB script which calculates both biased and unbiased ACF estimates
 of a signal and then use these ACF estimates to compute the corresponding
 correlogram in Eq. (15). Validate your code for different signals
e.g. WGN, noisy sinusoidal signals and ?ltered WGN. 

Explain how the spectral estimates based on (16)-(17) differ from one another? 
In particular, how does the correlogram corresponding to the unbiased ACF
 estimates behave for large lags (i.e. k close to N)? 

Does the unbiased ACF estimate result in negative values for the estimated PSD?
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 1;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
SAVE_FIGS = true;


%% Initialise Test Signals
close all % close current figures
fH = []; % clear the figure handle variable

% Test Signal Length
N = 2^10;

% Moving Average Filter - Transfer Function
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize); % numerator
a = 1;  % denominator

% Test Signal
testSignal.WGN = wgn(N, 1, 1); % white Gaussian Noise
testSignal.SINE = sin(linspace(0, N, N))' + random('Normal', 0, 1, [N, 1]); % random from: https://uk.mathworks.com/help/stats/prob.normaldistribution.random.html?s_tid=doc_ta
testSignal.WGN_F = filter(b, a, testSignal.WGN); % applying MA Filter
testSignal.types = {"WGN","SINE","WGN_F"};
testSignal.Names = {"WGN", "Noisy Sine", "Filtered WGN"};

for ii = 1:size(testSignal.Names,2)
    
    % Use Custom Func to find ACF & PSD
    [r_biased, r_unbiased, lags, PSD_biased, PSD_unbiased, fs] = ...
                            acfEstimator(testSignal.(testSignal.types{ii}));    
    fH{length(fH)+1} = figure;
        hold on;
        plot(lags, r_unbiased, 'DisplayName', 'unbiased');
        plot(lags, r_biased, 'DisplayName', 'biased');
        xlabel("Normalised Frequency, $\omega$");
        ylabel("PSD, $P(\omega)$");
        title(sprintf("\\textbf{%s}: Autocorrelation Function", testSignal.Names{ii}));
        legend('show');
        
    fH{length(fH)+1} = figure;
        hold on;
        plot(fs, PSD_unbiased, 'DisplayName', 'unbiased');
        plot(fs, PSD_biased, 'DisplayName', 'biased');
        xlabel("Lags, $k$");
        ylabel("ACF, $r(k)$");
        title(sprintf("\\textbf{%s}: Correlogram", testSignal.Names{ii}));
        legend('show');
end
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_3a_fig',num2str(ii,'%02i')],'pdf')
    end
end
