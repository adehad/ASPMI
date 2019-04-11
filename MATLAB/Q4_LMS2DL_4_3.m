%% Q4 From LMS to Deep Learning
% 4.3
%{
The activation can be generalised by scaling the activation function,
 that is, changing its amplitude (i.e. a·tanh). 
What range of values for a would be appropriate for the data in Figure 9? 
Pick a value of a in your suggested range and repeat the one-step ahead 
prediction with your new proposed activation function. 

Plot the zero-mean signal and the new prediction. 

Comment on the prediction and MSE & Rp in comparison with the standard LMS.
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
Delta = 1; % delta for prediction 

gamma = 0;
dPerc.activatorFunc = @tanh;
dPerc.bias = 0;
dPerc.ampl = max(y);

%% LMS
[U,~] = arima2diffEqns( y, lags(M)+Delta );  
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
        saveas(fH{ii},['figures', filesep,'q4_3_fig',num2str(ii,'%02i')],'pdf')
    end
end
