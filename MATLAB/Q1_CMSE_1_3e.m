%% Q1 Classical and Modern Spectrum Estimation
% 1.3e
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

% Test Signal Length
K = 50; % Total Length incl. Zero Pad

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
nSamples = 29;
numRealisations = 100;
rngSeed = 13; % initialise random seed
rng(rngSeed); % Set rng seed for consistency

sigSpaceDim = [1 2 3]; % Signal Space Dimensionality

clear PSE % to prevent problems with using a struct

for jj=1:length(sigSpaceDim)
    for ii = 1:numRealisations
        % Generate new sample space
        n = 0:nSamples;

        % Create the noisy signal
        testSignalTemp = testSignal(n); 
        testSignalTemp = [ testSignalTemp zeros(1, K-length(testSignalTemp)) ]; % zero pad

        % Auto correlation matrix estimate
        [~, R] = corrmtx(testSignalTemp, 14, 'mod'); % mod = modified
        % MUSIC
        structName = ['p', num2str(sigSpaceDim(jj))]; % create a name for the struct
        [ PSE.(structName)(ii, :) , Fs ] = ...
                                pmusic(R, sigSpaceDim(jj), [], 1, 'corr');
    end
end

%% Plot
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable

% Normalised Axis
fs = 0: 1/K : 1 -1/K;

% Plot Realisations and Mean
for jj=1:length(sigSpaceDim)
    fH{length(fH)+1} = figure;
        hold on;
        structName = ['p', num2str(sigSpaceDim(jj))]; % create a name for the struct
        plot(Fs, (PSE.(structName)), 'Color', COLORS(6, :), 'LineWidth', 0.5);
        plot(Fs, mean(PSE.(structName)), 'Color', COLORS(1, :));
        xlim([0.25 0.4])
        xlabel("Frequency ($\pi$ radians)");
        ylabel("PSD");
        title(sprintf("%d MUSIC Estimate Realisations and Mean \n $n=%d$, $p=%d$",numRealisations,nSamples,sigSpaceDim(jj)));

    fH{length(fH)+1} = figure;
        plot(Fs, std(PSE.(structName)), 'Color', COLORS(2, :));
        xlim([0.25 0.4])
        xlabel("Frequency ($\pi$ radians)");
        ylabel("PSD");
        title(sprintf("%d MUSIC Estimate Realisations Standard Deviation \n $n=%d$, $p=%d$",numRealisations,nSamples,sigSpaceDim(jj)));
end
    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_3e_fig',num2str(ii,'%02i')],'pdf')
    end
end
