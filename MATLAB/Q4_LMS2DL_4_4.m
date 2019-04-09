%% Q4 From LMS to Deep Learning
% 4.4
%{
 Now, consider the original time-series in Figure 9 which exhibits a non-zero mean.
 In order to account for the mean automatically, we can add a bias input to
 our model. 
In other words, the output of the model can be mathematically described by:
 phi(wT x + b) 
where phi( ) is the activation function, 
w are the weights of the model, x is the input data and b is the bias. 

Hint: the bias can be implemented by considering the ‘augmented’ input
 to the algorithm as [1, x]^T . 

Plot the original signal (non-zero mean) and the non-linear one-step ahead 
prediction with bias.
 Comment on the convergence of this algorithm.
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
[U,~] = arima2diffEqns( y, lags(M)+Delta );  
[y_pred,err,~] = dPerceptron(U, y', mu, gamma, dPerc.activatorFunc, dPerc);

msErr = meansqr(y-y_pred');
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
    title( sprintf('Dynamic Perceptron Output Variations \n $\\mu=$%.0e, $M=%i$, $\\phi$=\\texttt{%s}, a=%.2f, bias=%.2f',mu,M,func2str(dPerc.activatorFunc),dPerc.ampl,dPerc.bias) );
    xlabel('Time Index, $n$');
    ylabel('Magnitude');
    grid minor;
    legend('show', 'Location','best','NumColumns',2)
    xlim([0, 350])
    

    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q4_4_fig',num2str(ii,'%02i')],'pdf')
    end
end
