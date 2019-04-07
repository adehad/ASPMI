%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
% 3.1b
%{
Load the bivariate wind data of the wind speeds in the East-West direction, veast,
 and North-South direction, vnorth.
Form a complex-valued wind signal
v[n] = veast[n] + jvnorth[n] (38)
for the three wind regimes (low, medium, high). 

Use the scatter(x,y) function to plot a scatter diagram
 (the scatter diagram of the real and imaginary parts of a signal is
 also referred to as a circularity plot) for these three regimes. 

Comment on the circularity of the complex wind signal for low, medium 
and high wind speeds.

Configure the CLMS and ACLMS ?lters in a prediction setting (see Equation (16))
to perform a one-step ahead prediction of the complex wind data. 

Experiment with different ?lter lengths and comment on which algorithm (CLMS or
ACLMS) performs better for the different wind regimes.
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

% load wind data
load(['resources',filesep,'wind-dataset',filesep,'low-wind.mat']);
wind.low = complex(v_east, v_north);
load(['resources',filesep,'wind-dataset',filesep,'medium-wind.mat']);
wind.med = complex(v_east, v_north);
load(['resources',filesep,'wind-dataset',filesep,'high-wind.mat']);
wind.high = complex(v_east, v_north);

% number of samples
N = length(wind.low);
% sampling rate
Fs = 50e3; % sampled at 32e3, resampled to 50e3

windTypes = fieldnames(wind);

% Lecture 6, Slide 11:
circularity = @(data) abs( mean(data*data.') /mean(data*data') ); % NOTE: {.'} vs {'} operators
% I would have thought a better one would be to find major and minor axis
% by PCA and compare magnitude


% step-size - different for each wind type
mu = [0.1 0.01 0.001]; 
% model orders
M = 1:30;
lags = @(order) 0:order-1; % starts at 0 lag, hence -1

msErr = @(x) mean(abs(x).^2);

%% Circularity

rho = structfun(circularity,wind);

%% LMS Error
t_steady = 60; % steady state time index

for jj=1:length(windTypes)
    x = wind.(windTypes{jj});
    for ii = 1:length(M)
        % create differential equations
        [X, y] = arima2diffEqns(x, lags(M(ii)),1);

        % CLMS
        [~, err{1,ii}(:, :,jj), weights{1,ii}(:, :, jj)] = CLMS(X, x.', mu(jj));
        % ACLMS
        [~, err{2,ii}(:, :,jj), weights{2,ii}(:, :, jj)] = CLMS(X, x.', mu(jj), 'aug');

    end
end
    
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

windLabels = {'Low','Medium','High'};

% circularity investigation
for ii=1:length(windTypes)
    fH{length(fH)+1} = figure; hold on
        scatter(real(wind.(windTypes{ii})), imag(wind.(windTypes{ii})), 30, COLORS(ii,:), 'filled',  'DisplayName', windLabels{ii})
        title( sprintf('%s Regime \n Data Circularity $|\\rho|$=%.4f',windLabels{ii},rho(ii) ))
        xlabel("$v_{east}$");
        ylabel("$v_{north}$");
        legend('show','Location','best')
        grid minor
end

fH{length(fH)+1} = figure; hold on
    for ii=length(windTypes):-1:1
            scatter(real(wind.(windTypes{ii})), imag(wind.(windTypes{ii})), ...
                    20, COLORS(length(windTypes)-ii+1,:), 'filled', ...
                    'DisplayName', sprintf('%s: $|\\rho|$=%.3f',windLabels{ii},rho(ii)))
    end
    title('Comparison of the Wind Regimes')
    xLims = xlim;
    xlim([xLims(1), xLims(2)*1.5])
    xlabel("$v_{east}$");
    ylabel("$v_{north}$");
    legend('show','Location','best')
    grid minor


meanSqErr = cell2mat(cellfun(msErr,err,'UniformOutput',false));
% learning curves
for ii=1:length(windTypes)
    fH{length(fH)+1} = figure; hold on
        plot(M,pow2db(meanSqErr(1,:,ii)), 'DisplayName', 'CLMS');
        plot(M,pow2db(meanSqErr(2,:,ii)), 'DisplayName', 'ACLMS');
        title(sprintf('%s Regime \n Learning Curves of CLMS \\& ACLMS',windLabels{ii}))
        xlabel('Model Order')
        ylabel('MSE (dB)')
        legend('show','Location','best')
        grid minor
end
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_1b_fig',num2str(ii,'%02i')],'pdf')
    end
end
