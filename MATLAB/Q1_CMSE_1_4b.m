%% Q1 Classical and Modern Spectrum Estimation
% 1.4b
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


%% Initialise Auto Regressive Process

rngSeed = 13; % intialise random seed
rng(rngSeed); % enforce seed for consistency/reproducability

% Final Signal Length
N = 500;
Noffset = 500; % offset for calulation

% AR 'a' coefficients
a_AR = [2.76, -3.81, 2.65, -0.92];

% AR process
% randn provides the zero mean, 1 variance noise distribution
AR_Process = filter(1, [1 -a_AR], randn(N+Noffset,1) );
% Remove filter Transients by discarding first 500 samples
AR_Process = AR_Process(end-N:end); % end-N, enforces N long

% Ideal Frequency Response
[idealFreqResponse, w] = freqz(1, [1 -a_AR], N);

%% Error Estimating

% AR processes orders to try
modelOrder = 2:14;

% Pre-allocate for speed
  % error
    sigma_w = -ones(1, length(modelOrder)); % noise variance estimate
    error = zeros(1, length(modelOrder));   % error compared to ideal
  % log-likelihood
    logL = zeros(length(modelOrder), 1);    % ?

for ii = 1:length(modelOrder)
    
    % get coefficients with the modified covariance method: https://uk.mathworks.com/help/signal/autoregressive-and-moving-average-models.html
    [a_hat, sigma_w(1, ii)] = armcov(AR_Process, modelOrder(ii));
    
    [peaks(:, ii), ~] = freqz(sigma_w(1, ii)^(1/2), a_hat, N);
    
    error(:,ii) = mean(abs(peaks(:, ii)-idealFreqResponse).^2);
end
%%
%% Plot
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable

fH{1} = figure(1);   
    for ii = 1:length(modelOrder)
        plot(w/pi, pow2db(abs(idealFreqResponse).^2), "DisplayName", "ideal");
        hold on
        plot(w/pi, pow2db(abs(peaks(:, ii)).^2), "DisplayName", "model");
        title(sprintf("\\textbf{AR(%d)}: Spectral Estimation, $\\mathbf{N = %d}$", modelOrder(ii), N));
        xlabel("Normalised Frequency");
        ylabel("Power Spectral Density (dB)");
        grid on; grid minor;
        axis([0.15, 0.35, 0, 50]);
        legend("show");
        hold off
%         drawnow; pause; % comment out for interactive-ish mode
    end
    
fH{2} = figure;  
    plot(modelOrder, pow2db(sigma_w));
    title(sprintf("\\textbf{AR} Process: Noise Power against Model Order, $\\mathbf{N = %d}$", N));
    xlabel("Model Order, $p$");
    ylabel("Noise Power, $\sigma_{w}^{2}$ (dB)");
    grid on; grid minor;
    
fH{3} = figure;  
    plot(modelOrder, error);
    title(sprintf("\\textbf{AR} Process: MSE against Model Order, $\\mathbf{N = %d}$", N));
    xlabel("Model Order, $p$");
    ylabel("Mean Square Error $");
    grid on; grid minor;
%%
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
        title(sprintf("MUSIC Estimate: %d Realisations and Mean \n $n=%d$, $p=%d$",numRealisations,nSamples,sigSpaceDim(jj)));

    fH{length(fH)+1} = figure;
        plot(Fs, std(PSE.(structName)), 'Color', COLORS(2, :));
        xlim([0.25 0.4])
        xlabel("Frequency ($\pi$ radians)");
        ylabel("PSD");
        title(sprintf("MUSIC Estimate: %d Realisations Standard Deviation \n $n=%d$, $p=%d$",numRealisations,nSamples,sigSpaceDim(jj)));
end
    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_3e_fig',num2str(ii,'%02i')],'pdf')
    end
end
