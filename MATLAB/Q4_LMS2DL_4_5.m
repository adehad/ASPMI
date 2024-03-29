%% Q4 From LMS to Deep Learning
% 4.5
%{
Since we are only performing a singe weight update per time-step,
 it may take a lot of samples to converge to a reasonable prediction.
 
One solution to this issue is to pre-train the weights by over-?tting to a
 small number of samples. 
Start with w(0) = 0 and use 100 iterations (also called epochs) to
 fit the model to the first 20 samples to yield w_init.
 
Then use w_init as an initialisation to predict the entire time-series. 

Plot the original signal and the non-linear one-step ahead prediction,
 produced with bias and pre-trained weights.
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
y = y_raw;

% step-sizes
mu = 1e-5;
M = 4; % LMS order
lags = @(order) 0:order-1; % starts at 0 lag, hence -1
Delta = 1; % delta for prediction 

gamma = 0;
dPerc.activatorFunc = @tanh;
dPerc.bias = 1;
dPerc.ampl = max(y);

%% LMS
nSamples = 20;

[U,~] = arima2diffEqns( y(1:nSamples), lags(M)+Delta ); 
dPerc.W = zeros(size(U,1)+1, 1); % +1 as bias adds an extra dimension
% pre-train weights - overfit to first nSamples 
for ii=1:100
    [~,~,tempWeights] = dPerceptron(U, y(1:nSamples)', mu, gamma, dPerc.activatorFunc, dPerc);
    dPerc.W = tempWeights(:,end);
end
[U,~] = arima2diffEqns( y, lags(M)+Delta );  
dPerc.W = ones(size(U,1)+1,size(U,2)).*tempWeights(:,end);
[y_pred,err,~] = dPerceptron(U, y', mu, gamma, dPerc.activatorFunc, dPerc);
    
msErr(1) = meansqr(y-y_pred');
r_p(1) = mag2db(std(y_pred)/std(err));

fprintf('MSE: %.8f \t R_p: %.8f \n',msErr(1),r_p(1));

t_steady = 500;
msErr(2) = meansqr( y(t_steady:end)-y_pred(t_steady:end)' );
r_p(2) = mag2db( std( y_pred(t_steady:end) )/std( err(t_steady:end) ) );

fprintf('MSE: %.8f \t R_p: %.8f \n',msErr(2),r_p(2));
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% plot output
fH{length(fH)+1} = figure; hold on
    plot(y, 'DisplayName', 'True');
    plot(y_pred, 'LineStyle', ':', 'DisplayName', 'Prediction');
    title( sprintf('Dynamic Perceptron Output Variations \n $\\mu=$%.0e, $M=%i$, $\\phi$=\\texttt{%s}, a=%.2f, bias=%.2f',mu,M,func2str(dPerc.activatorFunc),dPerc.ampl,dPerc.bias) );
    xlabel('Time Index, $n$');
    ylabel('Magnitude');
    grid minor;
    legend('show', 'Location','best','NumColumns',2)
    xlim([0, 350])
    

    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q4_5_fig',num2str(ii,'%02i')],'pdf')
    end
end
