%% Q2 Adaptive Signal Processing
% 2.3b
%{
 Generate 1000 samples of s(n) and use your MATLAB program to estimate x(n)
for filter orders M = 5, 10, 15, 20.
Use values for Delta that range from the minimum value determined in part (a)
 to Delta = 25.
 What is the dependence (if any) between the delay Delta and the ' 
(i.e. how well the de-noising algorithm has performed)?

What is the effect of the filter order M on the MSPE? 

Given that the computational cost increases with M, what filter order
 would be a good pragmatic choice?
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
omega = 0.01/2;
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

clean_signal = @() ampl*sin(2*pi*omega*(1:N));


% step-sizes
mu = 0.01; % 0.05];

ALE.Delay = 0:25;

LMSorder = 5:5:20;
lags = @(order) 0:order-1; % starts at 0 lag, hence -1

for kk=1:length(LMSorder)
    % error
    err.diff{kk} = zeros(R, N, length(mu));
    % rows= realisations, cols = time index, pages = step size
    % MSE: steady-state
    err.mse{kk} = zeros(R, length(ALE.Delay));
    % MSPE: steady state
    err.mspe{kk} = zeros(1, length(ALE.Delay));
    
    for ii=1:length(ALE.Delay)
        x_hat{ii,kk} = zeros(R, N);
    end
end

    % rows= realisations, cols = time index, pages = coefficients
    % cell element = step size


%% Plot MSPE Error with Delay
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable
    
t_steady = 600; % steady state time index

omega_0 = [0.01/2, 0.01]; % correct frequency and double that

for omega=omega_0
    
    x = clean_signal();

    for kk=1:length(LMSorder)
        for ii=1:length(ALE.Delay)
            for jj=1:R
                s(jj,:) = noisy_signal(N,noise_mu,sigma_sq,b,a);

                % Convert process time series to differential eqn form
                [U,~] = arima2diffEqns( s(jj,:), lags(LMSorder(kk)), ALE.Delay(ii) );
                % LMS Error
                [x_pred,~,~] = LMS(U, s(jj,:), mu, 0);
                x_hat{ii,kk}(jj, :) = x_pred';

                err.diff{kk}(jj,:, ii) =  x - x_pred;

            end
            % MSE: steady-state
            err.mse{kk}(:,ii) = mean( err.diff{kk}(:, t_steady:end, ii).^2, 2); % mean over steady state
        end
    end
    % MSPE
    for kk=1:length(LMSorder)
        err.mspe{kk} = mean(err.mse{kk},1);
    end

    %% Plots
    % Plot MSPE variation with delay
    fH{length(fH)+1} = figure; hold on
        for kk=1:length(LMSorder)
            plot( ALE.Delay,  pow2db(err.mspe{kk}) , 'DisplayName', sprintf('$M$=%i',LMSorder(kk)));
        end

        title(sprintf('MSPE Error\n$\\mu=%.3f$, $\\omega_0=%.3f\\pi$',mu, omega));
        xlabel('Additional Delay, $\Delta$');
        ylabel('MSPE (dB)');
        grid minor;
        legend('show', 'NumColumns', 2, 'Orientation', 'horizontal')

end

%% Plot MSPE variation with model order
LMSorder = 0:20;
ALE.Delay = 3;

clear err
for kk=1:length(LMSorder)
    % error
    err.diff{kk} = zeros(R, N, length(mu));
    % rows= realisations, cols = time index, pages = step size
    % MSE: steady-state
    err.mse{kk} = zeros(R, length(ALE.Delay));
    % MSPE: steady state
    err.mspe{kk} = zeros(1, length(ALE.Delay));
    
    for ii=1:length(ALE.Delay)
        x_hat{ii,kk} = zeros(R, N);
    end
end

for kk=1:length(LMSorder)
    for ii=1:length(ALE.Delay)
        for jj=1:R
            s(jj,:) = noisy_signal(N,noise_mu,sigma_sq,b,a);

            % Convert process time series to differential eqn form
            [U,~] = arima2diffEqns( s(jj,:), lags(LMSorder(kk)), ALE.Delay(ii) );
            % LMS Error
            [x_pred,~,~] = LMS(U, s(jj,:), mu, 0);
            x_hat{ii,kk}(jj, :) = x_pred';

            err.diff{kk}(jj,:, ii) =  x - x_pred;

        end
        % MSE: steady-state
        err.mse{kk}(:,ii) = mean( err.diff{kk}(:, t_steady:end, ii).^2, 2); % mean over steady state
    end
end
% MSPE
for kk=1:length(LMSorder)
    err.mspe{kk} = mean(err.mse{kk},1);
end

%% Plots
% Plot MSPE variation with delay
fH{length(fH)+1} = figure; hold on
    plot( LMSorder,  pow2db(cell2mat(err.mspe)) );

    title(sprintf('MSPE Error\n$\\mu=%.3f$, $\\omega_0=%.3f\\pi$, $\\Delta$=%i',mu, omega, ALE.Delay));
    xlabel('Model Order, $M$');
    ylabel('MSPE (dB)');
    grid minor;


%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_3b_fig',num2str(ii,'%02i')],'pdf')
    end
end
