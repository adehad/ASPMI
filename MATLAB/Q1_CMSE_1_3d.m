%% Q1 Classical and Modern Spectrum Estimation
% 1.3d
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
% SAVE_FIGS = true;


%% Initialise Test Signals

% Test Signal Length
K = 512; % Total Length incl. Zero Pad
n = linspace(0,10,256); % Computed length
effectiveFs = 1/(2*n(2)); % frequency of linspace

% Amplitude parameters
A(1) = 1;
A(2) = 1;

% Frequency parameters
fHz(1) = 0.3;
fHz(2) = 0.32;

% Noise parameters
noiseV = 0.2; % noise variance
noiseMean = 0; % noise mean

% Noise Lambda Function
noiseTerm = @(n) noiseV/sqrt(2)*( randn(size(n))+1j*randn(size(n)) ) + noiseMean;

% Test Signal Lambda Function
testSignal = @(n) [ A(1)*exp(1j*2*pi*n*fHz(1)) + ...
                    A(2)*exp(1j*2*pi*n*fHz(2)) + ...
                                                noiseTerm(n) ]; 


                                    
%% Generate PSD
rngSeed = 13; % initialise random seed
rng(rngSeed); % Set rng seed for consistency
nSamples = [20 40 45 50];
for ii = 1:length(nSamples)
    
    % Generate new sample space
    n = 0:nSamples(ii);
    
    % Create the noisy signal
    testSignalTemp = testSignal(n); 
    testSignalTemp = [ testSignalTemp zeros(1, K-length(testSignalTemp)) ]; % zero pad
    
    % Find PSD
    PSD(ii,:) = fft( testSignalTemp ) / length(n);    

end

%% Plot
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable

% Normalised Axis
fs = 0: 1/K : 1 -1/K;

% Plot Realisations and Mean
fH{1} = figure;
    hold on;
    for ii=1:size(PSD,1)
        plotH(ii) = plot(fs, pow2db(abs(PSD(ii,:))));
        legendString(ii) = sprintf("$n=%d$", nSamples(ii));
    end
    yLims = ylim;
    plot([fHz(1), fHz(1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fHz(2), fHz(2)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    xlabel("Frequency ($\pi$ radians)");
    ylabel("PSD (dB)");
    title(strcat("PSD of Complex Exponentials"));
    legend(plotH,legendString);
    xlim([0.2 0.6])
    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_3d_fig',num2str(ii,'%02i')],'pdf')
    end
end
