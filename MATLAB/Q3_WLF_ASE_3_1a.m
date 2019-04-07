%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
% 3.1a
%{
Generate a first-order widely-linear-moving-average process, WLMA(1),
 driven by circular white Gaussian noise, x(n) [10]
y(n) = x(n) + b1x(n - 1) + b2x*(n - 1)      x = N (0, 1) (37)
where b1 = 1.5 + 1j and b2 = 2.5 - 0.5j. 

Write a MATLAB function for the ACLMS and implement both the
CLMS and ACLMS in the system identification setting to identify the WLMA model 
in (37). 

Plot the learning curve, 10log|e(n)|^2, for the ACLMS and CLMS.
 
Comment on the steady state error of the ACLMS and CLMS.

Hint: To obtain a smoother learning curve, plot the ensemble average of the 
learning curve from 100 independent simulations.
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

% number of samples
N = 1e3;
% realisations
R = 100;

% LMS parameters
mu = 0.1;
M = 2;
lags = @(order) 0:order-1; % starts at 0 lag, hence -1

% Noise WLMA model parameters
sigma_sq = 1; % consistency of the noise
b = [1.5+1i, 2.5-0.5i];
a = 1;

% Take in a signal (complex wgn) and apply the two weights to it
WLMA = @(a,b,signal,order) signal + conv( [ zeros(order, 1) ; b(1) ], ...
                                            signal(1:end-order) ) ...
                                  + conv( [ zeros(order, 1) ; b(2) ], ...
                                            conj(signal(1:end-order)) );
%% LMS Error
t_steady = 60; % steady state time index

for ii=1:R
    % guassian noise realisation
    x = wgn(N, 1,sigma_sq, 'complex');
    % create widely linear moving average 
    y = WLMA(a,b,x,1);

    % create differential equations
    [X, ~] = arima2diffEqns(x, lags(M));

    % CLMS
    [~, err{1}(:, ii), weights{1}(:, :, ii)] = CLMS(X, y', mu);
    % ACLMS
    [~, err{2}(:, ii), weights{2}(:, :, ii)] = CLMS(X, y', mu, 'aug');

end


%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% comparison on argand daigram of noise and actual signal
fH{length(fH)+1} = figure; hold on
    scatter(real(y), imag(y), 30, 'filled', 'DisplayName', 'WLMA(1)')
    scatter(real(x), imag(x), 30, 'filled', 'DisplayName', 'WGN')
    title('Comparison of Data Circularity')
    xlabel('Real $\Re$')
    ylabel('Imaginary $\Im$')
    legend('show','Location','best')
    grid minor

% learning curves
fH{length(fH)+1} = figure; hold on
    plot(pow2db(mean( abs(err{1}).^2, 2)), 'DisplayName', 'CLMS');
    plot(pow2db(mean( abs(err{2}).^2, 2)), 'DisplayName', 'ACLMS');
    title('Learning Curves of CLMS \& ACLMS')
    xlabel('Time index')
    ylabel('MSE (dB)')
    legend('show','Location','best')
    grid minor
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_1a_fig',num2str(ii,'%02i')],'pdf')
    end
end
