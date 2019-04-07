%% Q2 Adaptive Signal Processing
% 2.1f
%{
Implement the leaky LMS algorithm (with different values for gamma and mu)
 to estimate the AR coeffcients of the signal given in Part a). 

Why do the weights of the leaky LMS algorithm converge to 
incorrect values for the parameters a1 and a2?
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
a(1) = 0.1;     % 1 sample time delay coefficient
a(2) = 0.8;     % 2 sample time delay coefficient
sigma_sq = 0.25;
b = 0;          % AR numerator
p = length(a);  % model order

% AR process simulation
ar = arima('Constant', b, 'AR', a, 'Variance', sigma_sq);
x = simulate(ar, N, 'NumPaths', R);

% step-sizes
mu = [0.01 0.05];

% error
err.diff = zeros(R, N, length(mu));
    % rows= realisations, cols = time index, pages = step size
% weights
weights = {};
for ii=1:length(mu)
    weights{ii} = zeros(R, N, p);
end
    % rows= realisations, cols = time index, pages = coefficients
    % cell element = step size
% MSE: steady-state
err.mse = zeros(R, length(mu));
% EMSE: steady state
err.emse = zeros(R, length(mu));
% MisAdjustment: steady state
err.misAdj = zeros(R, length(mu));

%% LMS Error
t_steady = 400; % steady state time index
gamma = [0.1,0.5,1];
clear err % fix errors with re-using structs

for kk=1:length(gamma)
    for ii=1:length(mu)
        for jj=1:R
            % converyt ar process time series to differential eqn form
            [X2, y2] = arima2diffEqns(x(:, jj), [1:2]);
            % LMS Error
            [~, err.diff(jj, :, ii), tempWeights] = LMS(X2, x(:, jj)', mu(ii), gamma(kk));
            weights{ii,kk}(jj, :, :) = tempWeights';

        end
        % MSE: steady-state
        err.mse{kk}(:,ii) = mean( err.diff(:, t_steady:end, ii).^2, 2); % mean over steady state
    end
end

for kk=1:length(gamma)
    % EMSE
    err.emse{kk} = err.mse{kk} - sigma_sq;
    % MisAdjustment
    err.misAdj{kk} = err.emse{kk} / sigma_sq;
    
    fprintf('MisAdj = %.4f +/- %.5f \t [mu=%.2f, gamma=%.2f]\n', mean(err.misAdj{kk}(:,1)),std(err.misAdj{kk}(:,1)), mu(1), gamma(kk))
    fprintf('MisAdj = %.4f +/- %.5f \t [mu=%.2f, gamma=%.2f]\n', mean(err.misAdj{kk}(:,2)),std(err.misAdj{kk}(:,2)), mu(2), gamma(kk))
end




%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

whatMu = 1;
% whatGamma = 1;

for whatGamma=1:length(gamma)
    selWeightsA = weights{whatMu,whatGamma};
    selWeightsB = weights{whatMu+1,whatGamma};

    % Plot misAdjustment variation with time
    fH{length(fH)+1} = figure; hold on
        plot( err.misAdj{whatGamma}(:,1), 'DisplayName', sprintf('$\\mu=%.2f$',mu(1) ) );
        plot( err.misAdj{whatGamma}(:,2), 'DisplayName', sprintf('$\\mu=%.2f$',mu(2) ) );
        hold off

        title(sprintf('MisAdjustment Error \n $\\gamma=%.2f$', gamma(whatGamma) ));
        xlabel('Realisation');
        ylabel('MisAdjustment');
        grid minor;
        legend('show', 'Location', 'best', 'NumColumns', 2)

    % Plot Weight stabilisation over time
    fH{length(fH)+1} = figure;
        % mu = 0.05
        yyaxis left
        hold on
        plot( squeeze(mean(selWeightsA(:, :, 1))), 'DisplayName', '$\hat{a_{1}}$', 'Color', COLORS(1,:) );
        xLims{1} = xlim;
        plot( xlim, [a(1), a(1)], 'DisplayName', '${a_{1}}$', 'Color', COLORS(1,:), 'LineStyle', ':' );
        hold off

        yyaxis right
        hold on
        plot( squeeze(mean(selWeightsA(:, :, 2))), 'DisplayName', '$\hat{a_{2}}$', 'Color', COLORS(2,:) );
        plot( xlim, [a(2), a(2)], 'DisplayName', '${a_{2}}$', 'Color', COLORS(2,:), 'LineStyle', ':' );
        hold off

        yyaxis left
        title(sprintf('Weight Stabilisation with Time \n $\\mu=%.2f$ $\\gamma=%.2f$',mu(whatMu), gamma(whatGamma) ));
        xlabel('Time');
        ylabel('Weight Value', 'Color', 'k');
        grid minor;
        yyaxis right
        yLims = ylim;
        ylim([-inf, yLims(2)*1.2])
        legend('show','NumColumns',2,'Location','South')

    % Plot Weight stabilisation over time
    fH{length(fH)+1} = figure;
        % mu = 0.01
        yyaxis left
        hold on
        plot( squeeze(mean(selWeightsB(:, :, 1))), 'DisplayName', '$\hat{a_{1}}$', 'Color', COLORS(1,:) );
        xLims{1} = xlim;
        plot( xlim, [a(1), a(1)], 'DisplayName', '${a_{1}}$', 'Color', COLORS(1,:), 'LineStyle', ':' );
        hold off

        yyaxis right
        hold on
        plot( squeeze(mean(selWeightsB(:, :, 2))), 'DisplayName', '$\hat{a_{2}}$', 'Color', COLORS(2,:) );
        plot( xlim, [a(2), a(2)], 'DisplayName', '${a_{2}}$', 'Color', COLORS(2,:), 'LineStyle', ':' );
        hold off

        yyaxis left
        title(sprintf('Weight Stabilisation with Time \n $\\mu=%.2f$ $\\gamma=%.2f$',mu(whatMu+1), gamma(whatGamma) ));
        xlabel('Time');
        ylabel('Weight Value', 'Color', 'k');
        grid minor;
        yyaxis right
        yLims = ylim;
        ylim([-inf, yLims(2)*1.2])
        legend('show','NumColumns',2,'Location','South')

end

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_1f_fig',num2str(ii,'%02i')],'pdf')
    end
end
