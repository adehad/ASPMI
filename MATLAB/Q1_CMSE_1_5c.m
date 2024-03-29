%% Q1 Classical and Modern Spectrum Estimation
% 1.5c
%{
Plot the AR spectrum estimate for the RRI signals for the three trials. 
To find the optimal AR model order, experiment with your model order until
you observe a peak in the spectrum (approximately) corresponding to
the theoretical respiration rate. 

List the differences you observe between your estimated AR spectrum and the
periodogram estimate in Part a).
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
%{
% Only need to do once
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

% AR Process Model Orders
modelOrder = [2,4:6:24];

%% Periodograms
% Reset structs
PSD_standard = {};
AR_Coeff_Estimates = {};
PSD_AR_estimate = {};
sigma_w = {};
fAx = {};

for ii=1:size(RRI.x,2) % for all trials
    
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
    
    % Auto Regressive (AR) Process Estimate
     for jj=1:length(modelOrder)
        % get coefficients with the modified covariance method: https://uk.mathworks.com/help/signal/autoregressive-and-moving-average-models.html
        [AR_Coeff_Estimates{ii,jj}, sigma_w{ii,jj}] ...
                = armcov( tempRRI.*blackmanharris(N_orig), modelOrder(jj) );
     
        % PSD Estimate
        [PSD_AR_estimate{ii,jj}, fAx{ii}] ...
            = freqz(sqrt(sigma_w{ii,jj}), AR_Coeff_Estimates{ii,jj}, N_orig, RRI.fs);
     
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
    fH{length(fH)+1} = figure;
        plotH(1) = plot(fAx_standard, PSD_standard{ii},'DisplayName','standard');
        hold on
        for jj=1:length(modelOrder)
            plotH(1+jj) = plot( fAx{ii}, pow2db(abs(PSD_AR_estimate{ii,jj}).^2), 'DisplayName', sprintf('$p=%d$', modelOrder(jj)) );
        end
        
        % Harmonic plotter
        yLims = ylim;
        for jj=1:4 
             plot([fAx{ii}(find(fAx{ii}>=jj*trial_freqs(ii),1)), fAx{ii}(find(fAx{ii}>=jj*trial_freqs(ii),1))], [yLims(1), yLims(2)], 'LineWidth', 0.5, 'Color', [0 0 0 0.5], 'LineStyle','-.');
        end
        
        title(sprintf('Periodogram and AR Estimates: Trial %d', ii));
        xlabel('Frequency (Hz)');
        ylabel('PSD (dB/Hz)');
        grid on; grid minor;
        xlim([0 2]) % start at 0 to make one-sided , 0-2Hz is relevant range
%         ylim([-inf 0])
        legend(plotH,'NumColumns',2)
end

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_5c_fig',num2str(ii,'%02i')],'pdf')
    end
end
