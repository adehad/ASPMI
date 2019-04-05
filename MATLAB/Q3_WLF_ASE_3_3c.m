%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
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

% number of samples
N = 1.5e3;
% sampling rate
fs = 1.5e3; 

% noise
sigma_sq = 0.05;
eta = wgn(N, 1, pow2db(sigma_sq), 'complex');


% make frequencies
freq = zeros(1,N);  % pre-allocate
thresh = [5e2,1e3,1.5e3];        % thresholds
for n=1:length(freq)
    if n <= thresh(1)
        freq(n) = 100; 
    elseif n <= thresh(2)
        freq(n) = 100 + (n-500)/2;
    elseif n <= thresh(3)
        freq(n) = 100 + ( (n-1000)/25 )^2;
    else
        freq(n) = 0;
    end
end

% make phase - by integration
phi = cumtrapz(freq);

% freqz number of frequency bins on the unit circle
K = 2048;

% LMS
mu = 1; % step-size
gamma = [0 0.01 0.05 0.1 0.5]; % leakage factor


%% CLMS 
% NOTE: {.'} operator used to transpose, without complex conjugating
y = exp( (2*pi* phi/fs)* 1j ) + eta.';               % noisy signal
w = ( 0:(K-1) ).*(fs/K);                             % x axis  
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

% plot psd estimation        
for ii=1:length(clms_model.spectrum)
    fH{length(fH)+1} = figure; hold on
        surf(1:N, w, clms_model.spectrum{ii}, 'LineStyle','none','FaceColor','interp')
        view(2);
        title(sprintf("FM: DFT-AR(%i) Spectrogram $\\gamma=%.3f$", M, gamma(ii)));
        xlabel("Time Index");
        ylabel("Frequency (Hz)");
        grid on; grid minor;
        c = colorbar('eastoutside', 'TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
        c.Label.String = 'PSD (dB/Hz)';
        axis tight
end

    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_3c_fig',num2str(ii,'%02i')],'pdf')
    end
end
