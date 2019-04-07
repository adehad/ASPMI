%% Q1 Classical and Modern Spectrum Estimation
% 1.5a
%{
Apply the standard periodogram as well as the averaged periodogram with different
 window lengths (e.g. 50 s, 150 s ) to obtain the power spectral density 
of the RRI data. 

Plot the PSDs of the RRI data obtained from the three trials separately.
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


%% LOAD ECG DATA & CONVERT TO RRI
% Only need to do once
%{
load(['resources',filesep,'ECG_Data',filesep,'ECG_Split_110219.mat']);
% ECG: 3 cells for the 3 RRI recordings
% time_split: time axes for coresponding
% fs:   sampling  frequency

for ii=1:size(ECG,2)
    [RRI.x{ii}, RRI.fs] = ECG_to_RRI(ECG{ii},fs);
end
% Said yes to remove all anomolies

save(['resources',filesep,'ECG_Data',filesep,'RRI_Split_110219.mat'],'RRI');
%}
%% LOAD RRI DATA
load(['resources',filesep,'ECG_Data',filesep,'RRI_Split_110219.mat']);
% RRI.x{} - is a cell array of the RRI for each trial
% RRI.fs  - is the new sampling frequnecy 

% Recenter lambda function (removes mean bias)
recenter = @(x) (x - mean(x));

K = 2^12; % Zero-Padded Signal Length 

% Window Length (Seconds)
win{1} = 50;
win{2} = 150;
win{3} = 250;

%% Periodograms
% Reset structs
PSD_standard = {};
PSD_averaged = {};
PSD_avg_conf = {};
fAx = {};

for ii=1:size(RRI.x,2)
    
    % Load RRI Trial
    tempRRI = RRI.x{ii}';
    
    % Pre-process
    tempRRI = recenter(tempRRI);    % remove mean biases
    tempRRI = detrend(tempRRI);     % remove linear trends
    
    N_orig = length(tempRRI); % original data length
    
    % Standard Periodogram
    % DFT
    FFT_tempRRI = ( tempRRI .*blackmanharris(N_orig) ); % window (using blackman-harris)
    FFT_tempRRI = [ FFT_tempRRI' zeros(1, K-N_orig) ]; % zero padd the signal
    FFT_tempRRI = fftshift(fft(FFT_tempRRI)); % apply the fft
    PSD_standard{ii} = pow2db( abs(FFT_tempRRI).^2 /(N_orig *2*pi) )'; % convert to power spectrum (db)
    
    % Averaged Periodogram
    % Note: We use the Welch estimate and specify an overlap of 0 to achieve
    %       the Bartlett estimate
    overlap = 0; % samples of window overlap, if not specified defaults to 50% in pwelch
    for jj=1:size(win,2)
        [PSD_averaged{ii,jj}, fAx{ii}, PSD_avg_conf{ii,jj}] = ...
         pwelch(tempRRI, rectwin(win{jj}), overlap, K, RRI.fs, 'onesided', 'power','ConfidenceLevel',0.95);
    end
end

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable


fAx_standard = -1:(2/K):1-2/K; % normalised axis - NOTE: starts at -1 (currently not one sided periodogram)
fAx_standard = fAx_standard*RRI.fs/2;  % Hertz activate

trial_freqs = [0.3115,0.417,0.125];

% Standard Periodograms for all trials
for ii=1:size(RRI.x,2)
    fH{length(fH)+1} = figure; hold on
        plot(fAx_standard, PSD_standard{ii});

        % Harmonic plotter
        yLims = ylim;
        yLims(2) = 0;
        for jj=1:4 
             plot([fAx{ii}(find(fAx{ii}>=jj*trial_freqs(ii),1)), fAx{ii}(find(fAx{ii}>=jj*trial_freqs(ii),1))], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
        end
        
        title(sprintf('Standard Periodogram: Trial %d', ii));
        xlabel('Frequency (Hz)');
        ylabel('PSD (dB/Hz)');
        grid on; grid minor;
        xlim([0 2]) % start at 0 to make one-sided , 0-2Hz is relevant range
        ylim([-inf 0])
end


% Averaged Periodograms for all windows - for a given trial
for ii=1:size(RRI.x,2)
    fH{length(fH)+1} = figure;
        hold on
        for jj=1:size(win,2)
            % PSD
            plotH(jj) = plot(fAx{ii}, pow2db(PSD_averaged{ii,jj}),'Color',COLORS(jj,:));
            
            % Confidence Levels
%             plot(fAx{ii},pow2db(PSD_averaged_conf{ii,jj}),'Color',COLORS(ii,:),'LineWidth',0.5,'LineStyle','-.')
            
            legendString{jj} = sprintf('$W_L$ : %d', win{jj}); % unfortunately remakes it every loop
        end
        
        % Harmonic plotter
        yLims = ylim;
        for jj=1:4 
             plot([fAx{ii}(find(fAx{ii}>=jj*trial_freqs(ii),1)), fAx{ii}(find(fAx{ii}>=jj*trial_freqs(ii),1))], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
        end
        
        legend(plotH, legendString)
        title(sprintf('Bartlett Averaged Periodogram: Trial %d', ii));
        xlabel('Frequency (Hz)');
        ylabel('PSD (dB/Hz)');
        grid on; grid minor;
        xlim([0 2]) 
end

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_5a_fig',num2str(ii,'%02i')],'pdf')
    end
end
