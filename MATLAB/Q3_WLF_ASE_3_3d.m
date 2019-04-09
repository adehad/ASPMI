%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
% 3.3d
%{
Implement the DFT-CLMS for the EEG signal POz used in Part 1.2.
 To reduce computational burden, choose any segment POzof length 1200, 
e.g. POz(a:a+1200-1). 

Explain your observation about the time-frequency spectrum of the EEG signal.
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 3;
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

% Number of Samples
N = 3e3;

% length of DFT: 5 DFT samples per Hz
K = Fs*5;


% (noisy EEG) signal definition
startEl = 2e3;
y = s(startEl:startEl+N-1);
y = y - mean(y);            % remove signal dc offset, to reduce 0Hz component
% x axis  
w = ( 0:(K-1) ).*(Fs/K);                             


% LMS
mu = 1; % step-size
gamma = [0 0.001 0.005, 0.01]; % leakage factor


%% CLMS 
% NOTE: {.'} operator used to transpose, without complex conjugating
X = (1/K) *exp( 1j *(1:N)' *pi*( 0:(K-1) )/K ).';    % DFT 


for ii=1:length(gamma)
    
    % CLMS
    [~, ~,H] = CLMS(X, y, mu, gamma(ii));
    
    % power spectrum estimation
    H = abs(H).^2;

    % remove outliers
    medianH = 50 * median(median(H));
    outlier = H > medianH;
    H(outlier) = medianH;
    
    % store
    clms_model.spectrum{ii} = H;
    clms_model.w{ii} = w;
    clms_model.outlier(ii) = sum(outlier,'all');
end


%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

freqOfInterest = [13,26,50,100];

% plot psd estimation        
for ii=1:length(clms_model.spectrum)
    fH{length(fH)+1} = figure; hold on
        % plotting this as a 3D surf is unecessary, going to an image instead
        %{
        surf(1:N, w, clms_model.spectrum{ii}, 'LineStyle','none','FaceColor','interp')
        view(2);
        %}
        imagesc(1:N, w, clms_model.spectrum{ii}); axis xy;
        for jj=1:length(freqOfInterest)
            plot([1, N],[freqOfInterest(jj), freqOfInterest(jj)], 'r:');
        end
        title(sprintf("EEG: \\texttt{POz} DFT-CLMS Spectrogram $\\gamma=%.3f$", gamma(ii)));
        xlabel("Time Index");
        ylabel("Frequency (Hz)");
        grid on; grid minor;
        c = colorbar('eastoutside', 'TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
        c.Label.String = 'PSD (dB/Hz)';
        axis tight
        ylim([0 125])
        yticks(0:25:125);
end

    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_3d_fig',num2str(ii,'%02i')],'pdf')
    end
end
