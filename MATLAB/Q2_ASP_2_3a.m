%% Q2 Adaptive Signal Processing
% 2.3a
%{
Given the model for the noise in (34), what is the minimum value for the delay
 delta that may be used in the adaptive line enhancer in Fig 6 for a filter length 
of M > 1?
 
Hint: Start by expanding the Mean Square Error E{( s(n) - x^(n) )^2} to determine 
how the correlation in the noise affects the linear predictor.

Write a MATLAB program for the adaptive line enhancement using the LMS 
algorithm and justify your choice of minimum delta.
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


% step-sizes
mu = 0.01; % 0.05];

ALE.Delay = 0:20;

LMSorder = 5;
lags = @(order) 0:order-1; % starts at 0 lag, hence -1

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

for ii=1:length(ALE.Delay)
    for jj=1:R
        s(jj,:) = noisy_signal(N,noise_mu,sigma_sq,b,a);
        
        % Convert process time series to differential eqn form
        [U,~] = arima2diffEqns( s(jj,:), lags(LMSorder), ALE.Delay(ii) );
        % LMS Error
        [x_pred,~,~] = LMS(U, s(jj,:), mu, 0);
        x_hat{ii}(jj, :) = x_pred';
        
        err.diff(jj,:, ii) =  x - x_pred;
       
    end
    % MSE: steady-state
    err.mse(:,ii) = mean( err.diff(:, t_steady:end, ii).^2, 2); % mean over steady state
end

% MPSE
err.mpse = mean(err.mse,1);

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% Plot misAdjustment variation with time
fH{length(fH)+1} = figure;
    plot( ALE.Delay,  pow2db(err.mpse) );
%     hold on;
%     plot(diff(err.mpse));

    title('MPSE Error');
    xlabel('Additional Delay, $\Delta$');
    ylabel('MPSE (dB)');
    grid minor;


%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_3a_fig',num2str(ii,'%02i')],'pdf')
    end
end
