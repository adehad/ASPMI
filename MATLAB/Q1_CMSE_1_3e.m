%% Q1 Classical and Modern Spectrum Estimation
% 1.3e
%{
Use the following code to find the desired line spectra using the MUSIC method.
[X,R] = corrmtx(x,14,’mod’);
[S,F] = pmusic(R,2,[ ],1,’corr’);
plot(F,S,’linewidth’,2); set(gca,’xlim’,[0.25 0.40]);
grid on; xlabel(’Hz’); ylabel(’Pseudospectrum’);

Explain the operation of the first three lines in the code using 
the MATLAB documentation and the lecture notes. [10]
What is the meaning of the input arguments for the functions corrmtxand pmusic?
Does the spectrum estimated using the MUSIC algorithm provide more detailed
 information? 
State briefly the advantages and disadvantages of the periodogram and 
the MUSIC algorithms and comment on the bias and variance.
 
How accurate would a general spectrum estimate be when using MUSIC?
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
        xlabel('Frequency ($\pi$ radians)');
        ylabel('PSD (dB/Hz)');
        title(sprintf('%d MUSIC Estimate Realisations and Mean \n $n=%d$, $p=%d$',numRealisations,nSamples,sigSpaceDim(jj)));
        grid minor 
        
    fH{length(fH)+1} = figure;
        plot(Fs, std(PSE.(structName)), 'Color', COLORS(2, :));
        xlim([0.25 0.4])
        xlabel('Frequency ($\pi$ radians)');
        ylabel('PSD (dB/Hz)');
        title(sprintf('%d MUSIC Estimate Realisations Standard Deviation \n $n=%d$, $p=%d$',numRealisations,nSamples,sigSpaceDim(jj)));
        grid minor
end
    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_3e_fig',num2str(ii,'%02i')],'pdf')
    end
end
