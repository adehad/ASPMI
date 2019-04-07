%% Q2 Adaptive Signal Processing
% 2.3c
%{
Write a MATLAB program for the adaptive noise cancellation con?guration (see Fig 7).
 
Use the MSPE measure in (32) to compare the relative performance of the 
ANC and ALE in de-noising the noisy sinusoid in (33).
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

% number of samples
N = 1000;
% realisations
R = 100;

% Signal
% Sine Wave
ampl = 1;
omega = 0.01;
% noise + colouring filter
sigma_sq = 1;
noise_mu = 0;

% Filter Numerator
b(1) = 1;     % 0 sample time delay coefficient
b(2) = 0;     % 1 sample time delay coefficient
b(3) = 0.5;   % 2 sample time delay coefficient

a = 1;          % Filter Denominator
p = length(a);  % model order

noisy_signal = @(N,noise_mu,sigma_sq,b,a)  ampl*sin(2*pi*omega*(1:N)) ...
                    + filter( b, a, random('Normal', noise_mu, sigma_sq, 1, N) );
                
clean_signal = ampl*sin(2*pi*omega*(1:N));
                
noisyCorr.m = 0.8; % gradient for noise correlation
noisyCorr.c = +0.01; % offset for noise correlation

% step-sizes
mu = 0.01; % 0.05];
LMSorder = 6;
lags = 0:LMSorder-1; % starts at 0 lag, hence -1

ALE.Delay = 3;
ANC.Delay = 0; 


% error
err.diff = zeros(R, N, length(mu));
    % rows= realisations, cols = time index, pages = step size
% weights
weights = {};

for ii=1:length(ALE.Delay)
    x_hat{ii} = zeros(R, N);
end
    % rows= realisations, cols = time index, pages = coefficients
    % cell element = step size
% MSE: steady-state
err.mse = zeros(R, length(ALE.Delay));
% MisAdjustment: steady state
err.mpse = zeros(1, length(ALE.Delay));

%% LMS Error
t_steady = 300; % steady state time index
x = clean_signal;
filterTypes = {'ALE','ANC'};

for jj=1:R
    
    % Create Noisy signal
    s(jj,:) = noisy_signal(N,noise_mu,sigma_sq,b,a);
    eta_noise = s(jj,:) - x; % noisy - clean = just noise
    %% ALE
    % converyt ar process time series to differential eqn form
    [U,~] = arima2diffEqns( s(jj,:), lags, ALE.Delay );
    % LMS Error
    [x_pred,~,~] = LMS(U, s(jj,:), mu, 0);
    x_hat{1}(jj, :) = x_pred';

    err.diff(jj,:, 1) = x - x_pred;
    
    %% ANC
    % Create Correlated Noise
    eps_noise(jj,:) = noisyCorr.m*( eta_noise ) + noisyCorr.c;

    % converyt ar process time series to differential eqn form
    [U2,~] = arima2diffEqns( eps_noise(jj,:), lags, ANC.Delay );        
    % LMS Error
    [eta_pred,~,~] = LMS(U2, s(jj,:), mu, 0);
    x_hat{2}(jj, :) = s(jj,:) - eta_pred;

    err.diff(jj,:, 2) = x - x_hat{2}(jj, :);

end
% MSE: steady-state
err.mse = squeeze( mean( err.diff(:, t_steady:end, :).^2, 2) ); % mean over steady state

% MPSE
err.mpse = mean(err.mse,1);

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% plot mean predicted signal with time
fH{length(fH)+1} = figure;
    plot(  mean(x_hat{1}(:, :),1)  , 'DisplayName', 'ALE');
    hold on;
    plot(  mean(x_hat{2}(:, :),1) , 'DisplayName', 'ANC');
    plot( x, 'k:' , 'DisplayName', 'True');

    title('ALE vs ANC Mean');
    xlabel('Time Index');
    ylabel('Magnitude');
    grid minor;
    legend('show')
    
    
fH{length(fH)+1} = figure; hold on
    plotH(1,:) = plot(  x_hat{1}(:, :)' , 'Color',COLORS(1,:));
    plotH(2,:) = plot( x_hat{2}(:, :)' , 'Color',COLORS(2,:) );


    plotH(3,:) = plot( x , 'Color',COLORS(3,:) , 'LineStyle',':');

    title('ALE vs ANC Ensemble');
    xlabel('Time Index');
    ylabel('Magnitude');
    grid minor;
    legend(plotH(:,1), {'ALE','ANC','True'})


%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_3c_fig',num2str(ii,'%02i')],'pdf')
    end
end
