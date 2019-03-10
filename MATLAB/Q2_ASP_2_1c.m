%% Q2 Adaptive Signal Processing
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 2;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
% SAVE_FIGS = true;


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
ar = arima("Constant", b, "AR", a, "Variance", sigma_sq);
x = simulate(ar, N, "NumPaths", R);

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

for ii=1:length(mu)
    for jj=1:R
        % converyt ar process time series to differential eqn form
        [X2, y2] = arima2diffEqns(x(:, jj), [1:2]);
        % LMS Error
        [~, err.diff(jj, :, ii), tempWeights] = LMS(X2, x(:, jj)', mu(ii), 0);
        weights{ii}(jj, :, :) = tempWeights';
       
    end
    % MSE: steady-state
    err.mse(:,ii) = mean( err.diff(:, t_steady:end, ii).^2, 2); % mean over steady state
end

% EMSE
err.emse = err.mse - sigma_sq;
% MisAdjustment
err.misAdj = err.emse / sigma_sq;


fprintf('MisAdj = %.4f +/- %.5f \t [mu=%.2f]\n', mean(err.misAdj(:,1)),std(err.misAdj(:,1)), mu(1))
fprintf('MisAdj = %.4f +/- %.5f \t [mu=%.2f]\n', mean(err.misAdj(:,2)),std(err.misAdj(:,2)), mu(2))

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% Plot misAdjustment variation with time
fH{length(fH)+1} = figure;
    plot( err.misAdj(:,1), 'DisplayName', sprintf("$\\mu=%.2f$", mu(1)) );
    hold on
    plot( err.misAdj(:,2), 'DisplayName', sprintf("$\\mu=%.2f$", mu(2))  );
    hold off
    
    title("Mis-Adjustment Error");
    xlabel("Realisation");
    ylabel("Mis-Adjustment");
    grid minor;
    legend('show')

% Plot Weight stabilisation over time
fH{length(fH)+1} = figure;
    % mu = 0.05
    yyaxis left
    hold on
    plot( squeeze(mean(weights{1}(:, :, 1))), 'DisplayName', '$\hat{a_{1}}$', 'Color', COLORS(1,:) );
    xLims{1} = xlim;
    plot( xlim, [a(1), a(1)], 'DisplayName', '${a_{1}}$', 'Color', COLORS(1,:), 'LineStyle', ':' );
    hold off
    
    yyaxis right
    hold on
    plot( squeeze(mean(weights{1}(:, :, 2))), 'DisplayName', '$\hat{a_{2}}$', 'Color', COLORS(2,:) );
    plot( xlim, [a(2), a(2)], 'DisplayName', '${a_{2}}$', 'Color', COLORS(2,:), 'LineStyle', ':' );
    hold off
    
    yyaxis left
    title(sprintf("Weight Stabilisation with Time $\\mu=%.2f$",mu(1)));
    xlabel("Time");
    ylabel("Weight Value", 'Color', 'k');
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
    plot( squeeze(mean(weights{2}(:, :, 1))), 'DisplayName', '$\hat{a_{1}}$', 'Color', COLORS(1,:) );
    xLims{1} = xlim;
    plot( xlim, [a(1), a(1)], 'DisplayName', '${a_{1}}$', 'Color', COLORS(1,:), 'LineStyle', ':' );
    hold off
    
    yyaxis right
    hold on
    plot( squeeze(mean(weights{2}(:, :, 2))), 'DisplayName', '$\hat{a_{2}}$', 'Color', COLORS(2,:) );
    plot( xlim, [a(2), a(2)], 'DisplayName', '${a_{2}}$', 'Color', COLORS(2,:), 'LineStyle', ':' );
    hold off
    
    yyaxis left
    title(sprintf("Weight Stabilisation with Time $\\mu=%.2f$",mu(2)));
    xlabel("Time");
    ylabel("Weight Value", 'Color', 'k');
    grid minor;
    yyaxis right
    yLims = ylim;
    ylim([-inf, yLims(2)*1.2])
    legend('show','NumColumns',2,'Location','South')

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_1d_fig',num2str(ii,'%02i')],'pdf')
    end
end
