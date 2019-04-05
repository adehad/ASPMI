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
K = 1024;

% LMS
M = 1; % order
lags = @(order) 0:order-1; % starts at 0 lag, hence -1
mu = [0.001, 0.01, 0.05, 0.1, 0.2]; % step-size
Delta = 1;
%% CLMS 
% noisy signal
y = exp( (2*pi* phi/fs)* 1j ) + eta.';

% converyt ar process time series to differential eqn form
[X,~] = arima2diffEqns(y, lags(M), Delta );

for ii=1:length(mu)
    
    % CLMS
    [~, ~,clms_model.h(ii,:)] = CLMS(X, y, mu(ii));
    
    % spectrum pre-allocation 
    H = zeros(K, N);
    
    % power spectrum estimation
    for n = 1:N
        [h, w] = freqz(1, [1; -conj( clms_model.h(ii,n) )], K, fs);
        H(:, n) = abs(h).^2;
    end

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
        surf(1:N, clms_model.w{ii}, clms_model.spectrum{ii}, 'LineStyle','none','FaceColor','interp')
        view(2);
        title(sprintf("FM: CLMS-AR(%i) Spectrogram $\\mu=%.3f$", M, mu(ii)));
        xlabel("Time Index");
        ylabel("Frequency (Hz)");
        grid on; grid minor;
        c = colorbar('eastoutside', 'TickLabelInterpreter', 'latex', 'FontName', 'Palatino Linotype'); % oh why MATLAB ...
        c.Label.String = 'PSD (dB/Hz)';
        axis tight
end

% Percentage of the total number of estimates that were outliers
fH{length(fH)+1} = figure; hold on
        plot(mu, 100*clms_model.outlier./numel(clms_model.spectrum{ii}))
        title(sprintf("FM: CLMS-AR(%i) Outliers ", M));
        xlabel("$\mu$");
        ylabel("\% of Total Elements");
        grid minor;
        yLim = ylim;
        ylim([yLim(1), 1.1*yLim(2)])
    

    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_2b_fig',num2str(ii,'%02i')],'pdf')
    end
end
