%% Q2 Adaptive Signal Processing
% 2.3d
%{
Load a single-channel EEG Data5 (either Czor POz) from the file EEG_Data_Assignment2.mat.
 As seen in the previous assignments, there is a strong 50 Hz component introduced
 by the mains. Using the ANC con?guration, remove the 50 Hz component by
 generating a synthetic reference input composed of a sinusoid of 50 Hz corrupted
by white Gaussian noise. 
(Hint: Experiment with different step-sizes, mu, and filter lengths, M, to 
suppress the 50 Hz component without affecting the other frequency components).
 
Plot the spectrum for the corrupted and denoised EEG data using the
 spectrogram function with appropriate window length and overlap.
 What do you observe?
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 2;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
SAVE_FIGS = true;


%% LMS Parameter
% Set rng seed for consistency
rSeed = 13;
rng(rSeed)

% Load EEG Data
eeg = load(['resources',filesep,'EEG_Data',filesep,'EEG_Data_Assignment2.mat']);

% sampling frequency
Fs = eeg.fs;
% signal
s = eeg.POz';

% number of samples
N = length(s);

% Noise
% Sine Wave
ampl = 1;
omega = 50/Fs; % 50Hz
% noise 
sigma_sq = 0.1; % consistency of the noise
noise_mu = 0;

noisy_signal = @(N,noise_mu,sigma_sq)  ampl*sin(2*pi*omega*(1:N)) ...
                    + random('Normal', noise_mu, sigma_sq, 1, N) ;
                
clean_signal = ampl*sin(2*pi*omega*(1:N));
                
noisyCorr.m = 0.8; % gradient for noise correlation
noisyCorr.c = +0.01; % offset for noise correlation

% step-sizes
mu = [0.001, 0.005, 0.01, 0.05, 0.1];
M = [1,5,10,15,25]; % LMS order
lags = @(order) 0:order-1; % starts at 0 lag, hence -1

% Spectogram
spect.winSize = 2^12;   % window size
spect.win = blackman(spect.winSize); % actual window
spect.overlap = 0.8;        % percentage window overlap
spect.nOverlap = floor(spect.winSize*spect.overlap);
spect.nFFT = 10*Fs;          % number of fft samples, e.g. 10 per sample
%% LMS Error
t_steady = 60*Fs; % steady state time index

% Reference Signal (epsilon)
eps_noise = noisy_signal(N,noise_mu,sigma_sq);

% Original
[Pss, ws] = periodogram(s, blackman(N), spect.nFFT, Fs,'power', 'onesided');

for jj=1:length(M)
    %% ANC
    % converyt ar process time series to differential eqn form
    [U,~] = arima2diffEqns( eps_noise, lags(M(jj)), 1 );  
%     d = [0 s(1:end-1)];
    d = s;
    for ii=1:length(mu)
        % LMS Error
        [eps_pred,~,~] = LMS(U, d, mu(ii), 0);
        x_hat(ii, :, jj) = d - eps_pred;
        
        % Periodogram of Cleaned Signal
        [Pxx(ii,:,jj), ~] = periodogram( x_hat(ii, :, jj), blackman(N), spect.nFFT, Fs,'power','onesided');
        
        try
        err.db(ii,:,jj) = pow2db(Pss)'-pow2db(squeeze(Pxx(ii,:,jj)));
        catch
            warning('some bad happened here')
        end
        
    end
end

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% plot original spectogram with others
fH{length(fH)+1} = figure;
    spectrogram(s, spect.win, spect.nOverlap, spect.nFFT, Fs, 'xaxis');
    xticks(0:25:125);
    xlim([0 125]);
    title('Original Spectogram')
    c = colorbar('TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
    c.Label.String = 'PSD (dB/Hz)';

% Print all periodograms
%{
for jj=1:length(M)
    for ii=1:length(mu)
        fH{length(fH)+1} = figure;
%             spectrogram(squeeze(x_hat(ii,t_steady:end,jj)), spect.win, spect.nOverlap, spect.nFFT, Fs, 'yaxis');
            spectrogram(squeeze(x_hat(ii,:,jj)), spect.win, spect.nOverlap, spect.nFFT, Fs, 'xaxis');
            xticks(0:25:125);
            xlim([0 125]);
            title(sprintf('$\\mu$=%.4f, M(%d)',mu(ii),M(jj)))
            c = colorbar('TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
            c.Label.String = 'PSD (dB/Hz)';
    end
end
%}

% Selected Mu and M to plot
whatMu = [1,5];
whatM = [3];   
   
% plot all selected
for whatM=whatM
    for whatMu=whatMu
        selPxx = squeeze( Pxx((whatMu), :, (whatM)) );
        fH{length(fH)+1} = figure; hold on
            plot(ws, pow2db(Pss), 'DisplayName', '$P_{Oz}$');
            plot(ws, pow2db(selPxx), 'DisplayName', 'ANC');
            xticks(0:25:125);
            xlim([30 125]);
            title(sprintf('$P_{Oz}$ vs ANC Periodogram \n $\\mu$=%.4f, M(%d)',mu(whatMu),M(whatM)))
            xlabel('Frequency (Hz)');
            ylabel('PSD (dB/Hz)');
            grid minor
            legend('show', 'Location', 'northeast');
            
        fH{length(fH)+1} = figure;          
            spectrogram(squeeze(x_hat(whatMu,:,whatM)), spect.win, spect.nOverlap, spect.nFFT, Fs, 'xaxis');
            xticks(0:25:125);
            xlim([0 125]);
            title(sprintf('ANC Filtered Periodogram \n $\\mu$=%.4f, M(%d)',mu(whatMu),M(whatM)))
            c = colorbar('TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
            c.Label.String = 'PSD (dB/Hz)';
    end
end      


% Error Calculation
interestRange = ws<=50+2 & ws>=50-2; % 2Hz on either side
err.mse_db      = squeeze(mean(err.db(:,interestRange,:).^2,2)); % 50Hz error - i.e. attenuation
err.mse_db_no50 = squeeze(mean(err.db(:,interestRange,:).^2,2)); % error for non-50db
% fprintf('Mean Error (excl. 50Hz): %.3f dB \n', err.mse_db() )

% old custom colormap when was using image() instead of imagesc()
%     err.mse_max = max(err.mse_db,[],'all','omitnan');
%     err.mse_min = min(err.mse_db,[],'all','omitnan');
%     cMap = parula( round(err.mse_max-err.mse_min) );
%     cMap = [ones(floor(err.mse_min),3).*cMap(1,:) ; cMap];

% Plot hyperparmeter heatmap
fH{length(fH)+1} = figure; hold on
    subplot(1,2,1)
        im = imagesc(mu, M, sqrt(err.mse_db));
        ax = gca; ax.YDir = 'normal'; axis tight;
        title('50Hz RMSE Error (dB)');
        colorbar('eastoutside','TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
        xLab = xlabel('$\mu$'); set(xLab, 'Units', 'Normalized', 'Position', [1.4405,-0.05325,0]);
        ylabel('M');
        xticks(mu); xticks([min(mu):((max(mu)-min(mu))/(length(mu)-1)):max(mu)-min(mu),max(mu)]); xticklabels(num2str(mu')); xtickangle(45)
        yticks([min(M):(max(M)-min(M))/(length(M)-1):max(M)-min(M),max(M)]); yticklabels(num2str(M'));

    subplot(1,2,2)
        imagesc(mu, M, sqrt(err.mse_db_no50) );
        colorbar('eastoutside','TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
        ax = gca; ax.YDir = 'normal'; axis tight;
        title('Non-50Hz RMSE Error (dB)')
%         xlabel('$\mu$');
        ylabel('M');
        xticks(mu); xticks([min(mu):((max(mu)-min(mu))/(length(mu)-1)):max(mu)-min(mu),max(mu)]); xticklabels(num2str(mu')); xtickangle(45)
        yticks(M); yticks([min(M):(max(M)-min(M))/(length(M)-1):max(M)-1,max(M)]); yticklabels(num2str(M'));

%     [textX, textY] = meshgrid(mu,M);
%     text(textX(:), ...
%          textY(:), ...
%          string(num2str(err.mse_db(:))))
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_3d_fig',num2str(ii,'%02i')],'pdf')
    end
end
