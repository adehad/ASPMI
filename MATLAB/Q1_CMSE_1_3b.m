%% Q1 Classical and Modern Spectrum Estimation
% 1.3b
%{
Use your code from the previous section (only the biased ACF estimator) 
to generate the PSD estimate of several realisations of a random process
 and plot them as in Fig. 1.
 
Generate different signals composed of sinusoids corrupted by noise and
 elaborate on how disperse are the different realisation of the spectral estimate.
 
Hint: use the
fftand fftshiftcommands in MATLAB.
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
rSeed = 13; % initialise random seed
rng(rSeed)

numRealisations = 40; % Number of Realisations
rngList = randi(2^16,numRealisations,1); % define random number list for consistent plots
noisePower = 3; % Noise Power
PSD_biased = zeros(numRealisations, K*2-1); % Pre-allocate for speed

for ii = 1:numRealisations
    % Generate a new noise term (each loop)
    noiseTerm = wgn( length(testSignal),1, noisePower, [], rngList(ii) );
    
    % Create the noisy signal
    testSignalTemp = testSignal + noiseTerm'; % transpose noise so same dimensions
    
    % Find PSD
    [~,~,~, PSD_biased(ii,:), ~, fs] = acfEstimator(testSignalTemp); 
                               
end

%% Plot
close all % close current figures
fH = []; % clear the figure handle variable

% Plot Realisations and Mean
fH{1} = figure; hold on;
    plot(fs*effectiveFs, real(PSD_biased)', 'Color', COLORS(6, :), 'LineWidth', 0.5);
    plot(fs*effectiveFs, mean(real(PSD_biased),1)', 'Color', COLORS(1, :));
    
    yLims = ylim;
    plot([fs(find(fs>=fHz(1),1)-1), fs(find(fs>=fHz(1),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fs(find(fs>=fHz(2),1)-1), fs(find(fs>=fHz(2),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fs(find(fs>=fHz(3),1)-1), fs(find(fs>=fHz(3),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fs(find(fs>=fHz(4),1)-1), fs(find(fs>=fHz(4),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');


    xlabel('Frequency ($\pi$ radians)');
    ylabel('PSD, $P(\omega)$');
    title(strcat(num2str(numRealisations), ' PSD Realisations \& Mean'));
%     legend('show');
    xlim([0 2])
    grid minor;
    
    
% Plot Standard Deviations    
fH{2} = figure; hold on;
    plot(fs*effectiveFs, std(real(PSD_biased),1)', 'Color', COLORS(2, :));
    
    yLims = ylim;
    plot([fs(find(fs>=fHz(1),1)-1), fs(find(fs>=fHz(1),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fs(find(fs>=fHz(2),1)-1), fs(find(fs>=fHz(2),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fs(find(fs>=fHz(3),1)-1), fs(find(fs>=fHz(3),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
    plot([fs(find(fs>=fHz(4),1)-1), fs(find(fs>=fHz(4),1)-1)], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');

    xlabel('Frequency ($\pi$ radians)');
    ylabel('PSD, $P(\omega)$');
    title(strcat(num2str(numRealisations), ' PSD Realisations Standard Deviation'));
%     legend('show');
    xlim([0 2])
    grid minor
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_3b_fig',num2str(ii,'%02i')],'pdf')
    end
end
