%% Q4 From LMS to Deep Learning
% 4.1
%{
Load the time-series in Figure 9 from the ‘time-series.mat’.
 Remove the mean from the time-series and use the LMS algorithm 
(with a learning rate of mu = 1e-5) from Section 3 to predict y[n] from
y[n - 1], y[n - 2], y[n - 3] and y[n - 4] 
(i.e. assume the data is generated using an AR(4) process). 

Plot the zero-mean version of the signal, y, against its one-step ahead prediction. 
Describe the behaviour of the output signal,
and calculate the mean-square error (MSE) 
and prediction gain, Rp = 10 log10(sigma_y ^2 / sigma_e ^2)
sigma_y is the variance of the LMS output
sigma_e is the variance of the prediction error
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 4;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
SAVE_FIGS = true;


%% LMS Parameter
% Set rng seed for consistency
rSeed = 13;
rng(rSeed)

% Load time-series
y_raw = load(['resources',filesep,'time-series.mat'],'y');
y_raw = y_raw.y;

% remove mean
y = y_raw - mean(y_raw);

% step-sizes
mu = 1e-5;
M = 4; % LMS order
lags = @(order) 0:order-1; % starts at 0 lag, hence -1
gamma = 0;

%% LMS
[U,~] = arima2diffEqns( y, lags(M) );  
[y_pred,err,~] = dPerceptron(U, y', mu, gamma);

msErr = meansqr(y-y_pred);
r_p = mag2db(std(y)/std(err));

fprintf('MSE: %.8f \t R_p: %.8f \n',msErr,r_p);
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% plot output
fH{length(fH)+1} = figure; hold on
    plot(y, 'DisplayName', 'True');
    plot(y_pred, 'LineStyle', ':', 'DisplayName', 'Prediction');
    title(sprintf('LMS Output Variations \n $\\mu=$%e, $M=%i$ ',mu,M));
    xlabel('Time Index, $n$');
    ylabel('Magnitude');
    grid minor;
    legend('show', 'Location','best','NumColumns',2)
    xlim([0, 350])

    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q4_1_fig',num2str(ii,'%02i')],'pdf')
    end
end
