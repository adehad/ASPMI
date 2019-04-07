%% Q2 Adaptive Signal Processing
% 2.1b
%{
Implement an LMS adaptive predictor using N = 1000 samples of x(n) as in a)
 and plot the squared prediction error e2(n) dB i.e. 10 log(e2(n))
 along time using step sizes µ = 0.05 and µ = 0.01.
 
Repeat this experiment for 100 different realizations of x(n)
 and plot the learning curve by averaging the plots of 10 log(e2(n)).
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

% AR model parameters
a(1) = 0.1;
a(2) = 0.8;
sigma_sq = 0.25;
b = 0;
p = length(a);

% AR process simulation
ar = arima('Constant', b, 'AR', a, 'Variance', sigma_sq);
x = simulate(ar, N, 'NumPaths', R);

% step-sizes
mu = [0.01 0.05];

% error
err = zeros(R, N, length(mu));
% rows= realisations, cols = time index, pages = step size
%% LMS Error
for ii=1:length(mu)
    for jj=1:R
        % convert AR process time series to differential eqn form
        [X, y] = arima2diffEqns(x(:, jj), [1:2]);
        % LMS Error
        [~, err(jj, :, ii), ~] = LMS(X, x(:, jj)', mu(ii), 0);
    end
end

% squared error
err = err.^2; % element wise square

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable


% Plot Realisation 12
fH{length(fH)+1} = figure; hold on
    plot(pow2db( err(12,:,1) ), 'DisplayName', sprintf('$\\mu=%.2f$', mu(1)) );
    plot(pow2db( err(12,:,2) ), 'DisplayName', sprintf('$\\mu=%.2f$', mu(2))  );
    hold off
    
    title('Squared Prediction Error, Realisation 12');
    xlabel('Time Increment');
    ylabel('Squared Error (dB)');
    grid minor;
    legend('show')

% Mean of log Square Error Trace
fH{length(fH)+1} = figure; hold on
    plot(pow2db( mean(err(:,:,1)) ), 'DisplayName', sprintf('$\\mu=%.2f$', mu(1)) );
    plot(mean( pow2db(err(:,:,2)) ), 'DisplayName', sprintf('$\\mu=%.2f$', mu(2)) );
    hold off
    
    title('Mean of Log of Squared Prediction Error');
    xlabel('Time Increment');
    ylabel('Mean Squared Error (dB)');
    grid minor;    
    legend('show')

% Log of Mean Trace
fH{length(fH)+1} = figure; hold on
    plot(pow2db( mean(err(:,:,1)) ), 'DisplayName', sprintf('$\\mu=%.2f$', mu(1)) );
    plot(pow2db( mean(err(:,:,2)) ), 'DisplayName', sprintf('$\\mu=%.2f$', mu(2)) );
    hold off
    
    title('Log of Mean Squared Prediction Error');
    xlabel('Time Increment');
    ylabel('Mean Squared Error (dB)');
    grid minor;    
    legend('show')
    
    % intuitively error high as when we not in steady state at the weights
    % are sorting themselves out

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_1b_fig',num2str(ii,'%02i')],'pdf')
    end
end
