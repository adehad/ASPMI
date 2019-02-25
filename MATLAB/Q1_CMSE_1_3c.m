%% Q1 Classical and Modern Spectrum Estimation
% 1.3c
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
A(1) = 2;
A(2) = 1.75;
A(3) = 0.85;
A(4) = 1.2;
% Frequency parameters
fHz(1) = 0.4;
fHz(2) = 0.6;
fHz(3) = 0.85;
fHz(4) = 0.95;

% Test Signal - (NOTE: Noise added during realisations)
testSignal = [ A(1)*sin(2*pi*n*fHz(1)) + ...
               A(2)*sin(2*pi*n*fHz(2)) + ...
               A(3)*sin(2*pi*n*fHz(3)) + ...
               A(4)*sin(2*pi*n*fHz(4)) ...
                                        zeros(1, K-length(n))]; % Zero Pad

%% Generate PSD
rngSeed = 13; % initialise random seed
numRealisations = 40; % Number of Realisations
rngList = randi(2^16,numRealisations,1); % define random number list for consistent plots
noisePower = 3; % Noise Power
PSD_biased = zeros(numRealisations, K*2-1); % Pre-allocate for speed
for ii = 1:numRealisations
    
    % Generate a new noise term (each loop_
    noiseTerm = wgn( length(testSignal),1, noisePower, [], rngList(ii) );
    
    % Create the noisy signal
    testSignalTemp = testSignal + noiseTerm'; % transpose noise so same dimensions
    
    % Use Custom Func to find PSD
    [~,~,~, PSD_biased(ii,:), ~, fs] = ...
                            acfEstimator(testSignalTemp);    

end

%% Plot
close all % close current figures
fH = []; % clear the figure handle variable

% Plot Realisations and Mean
fH{1} = figure;
    hold on;
    plot(fs*effectiveFs, pow2db(real(PSD_biased))', 'Color', COLORS(6, :), 'LineWidth', 0.5);
    plot(fs*effectiveFs, pow2db(mean(real(PSD_biased),1))', 'Color', COLORS(1, :));
    xlabel("Frequency ($\pi$ radians)");
    ylabel("PSD (dB)");
    title(strcat("PSD ", num2str(numRealisations), " Realisations \& Mean"));
%     legend('show');
    xlim([0 2])
    
% Plot Standard Deviations    
fH{2} = figure;
    hold on;
    plot(fs*effectiveFs, std(pow2db(real(PSD_biased)),1)', 'Color', COLORS(2, :));
    xlabel("Frequency ($\pi$ radians)");
    ylabel("PSD (dB)");
    title(strcat("PSD ", num2str(numRealisations), " Realisations Standard Deviation"));
%     legend('show');
    xlim([0 2])
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        savefig(fH{ii},['figures', filesep,'q1_3c_fig',num2str(ii,'%02i')])
    end
end
