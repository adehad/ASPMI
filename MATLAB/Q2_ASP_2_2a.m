%% Q2 Adaptive Signal Processing
% 2.2a
%{
 Implement the three GASS algorithms above, and compare their performance
 when identifying a real-valued MA(1) system
x(n) = 0.9eta(n - 1) + eta(n) eta = N (0, 0.5). (27)

Explain the performance advantages/disadvantages against a standard LMS algorithm 
(with a fixed step size, mu = 0.01 and mu = 0.1).
 
Plot the weight error curves (w˜(n) = wo - w(n) where wo = 0.9) for each algorithm
 and comment on their convergence speed and steady state error.
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

% VSS parameters
rho = 5e-3;
mu_0 = 0.1;
gamma = 0;
alpha = 0.8;

% GASS Types
gassType = {'Benveniste', 'Ang_Farhang', 'Matthews_Xie'};

% Noise model parameters
noise_mu = 0;
sigma_sq = 0.5;

% Moving Average (1) parameters
a(1) = 1.0;     % 0 sample time delay coefficient
a(2) = 0.9;     % 1 sample time delay coefficient
b = 1;          % numerator
p = length(a);  % model order

% step-sizes
mu = [0.01 0.05 0.1];

% error
err.diff = {}; % zeros(R, N, length(mu));
    % rows= realisations, cols = time index, pages = step size
    
% weights
weights = {};

for jj=1:length(gassType)+1 % +1 for normal LMS case
    for ii=1:length(mu)
        err.diff{ii,jj} = zeros(R, N);
        weights{ii,jj} = zeros(R, N, p);
        err.weight{ii,jj} = zeros(length(a),N);
    end
    
end
    % rows= realisations, cols = time index, pages = coefficients
    % cell rows = step size, columns = GASS type
% MSE: steady-state
err.mse = zeros(R, length(mu));
% EMSE: steady state
err.emse = zeros(R, length(mu));
% MisAdjustment: steady state
err.misAdj = zeros(R, length(mu));



%% LMS vs NLMS (w/ VSS) Error
t_steady = 400; % steady state time index
% gamma = [0.1,0.5,1];
gamma = 0;

for jj=1:R
    % Create noise realisation
    noise_eta = random('Normal', noise_mu, sigma_sq, 1, N);
    
    % Moving Average Filter
    x = filter(a, b, noise_eta);
    
    [X,~] = arima2diffEqns(noise_eta, [0,1]);
    
    for ii=1:length(mu)
        % LMS Error - NOTE: uses end element of cell arrays
        [~, err.diff{ii,end}(jj, :), tempWeights] = LMS(X, x, mu(ii), gamma);
        weights{ii,end}(jj, :, :) = tempWeights';
        
        for kk=1:length(gassType) 
           % VSS
            [~, err.diff{ii,kk}(jj, :), tempWeights]  = LMS_VSS(X, x,  mu(ii), gamma, rho, string(gassType{kk}), alpha);
            weights{ii,kk}(jj, :, :) = tempWeights';

            % MSE: steady-state
%             err.mse(:,ii) = mean( err.diff(:, t_steady:end, ii).^2, 2); % mean over steady state
            
        end
    end
end

for jj=1:length(a)
    for ii=1:length(mu)
        for kk=1:length(gassType)+1 % +1 for LMS
            err.weight{ii,kk}(jj,:) = a(jj) - squeeze(mean(weights{ii,kk}(:, :, jj)));
%             err.sqrE{ii,kk}() = err.diff{ii,kk}(jj, :, ii);
        end
    end
end

% EMSE
% err.emse = err.mse - sigma_sq;
% MisAdjustment
% err.misAdj = err.emse / sigma_sq;


%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

gassTypeLabels =  {'Benveniste', 'Ang_Farhang', 'Matthews_Xie'};
gassTypeLabels{length(gassTypeLabels) + 1} = 'LMS';
gassTypeLabels = replace(gassTypeLabels,'_',' \& ');

whatMu = 3;
whatGass = 2;


selError = err.weight{whatMu,whatGass};
selWeightsA = weights{whatMu,whatGass};

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
    title(sprintf('Weight Stabilisation with Time \n %s : $\\mu=%.2f$',gassTypeLabels{whatGass}, mu(whatMu) ));
    xlabel('Time');
    ylabel('Weight Value', 'Color', 'k');
    grid minor;
    yLims = ylim;
    ylim([yLims(1), yLims(2)*1.2])
    yyaxis right
    ylim([yLims(1), yLims(2)*1.2])
    legend('show','NumColumns',2,'Location','South')

% Plot Error stabilisation over time
fH{length(fH)+1} = figure;
    % mu = 0.05
    yyaxis left
    hold on
    plot( selError(1, :), 'DisplayName', '$\hat{a_{1}}$', 'Color', COLORS(1,:) );
    xLims{1} = xlim;
%     plot( xlim, [a(1), a(1)], 'DisplayName', '${a_{1}}$', 'Color', COLORS(1,:), 'LineStyle', ':' );
    yLims = ylim;
    ylim([-inf, yLims(2)+0.1])
    hold off

    yyaxis right
    hold on
    plot( selError(2,:), 'DisplayName', '$\hat{a_{2}}$', 'Color', COLORS(2,:) );
%     plot( xlim, [a(2), a(2)], 'DisplayName', '${a_{2}}$', 'Color', COLORS(2,:), 'LineStyle', ':' );
    hold off
    ylim([-inf, yLims(2)+0.1])

    yyaxis left
    title(sprintf('Weight Error with Time \n %s : $\\mu=%.2f$',gassTypeLabels{whatGass}, mu(whatMu) ));
    xlabel('Time');
    ylabel('Weight Error', 'Color', 'k');
    grid minor;
%     yyaxis right

    legend('show','NumColumns',2,'Location','South')

    %%
% Plot Weight Error Curve over time - For all
lineStyleList = {'-',':','--'};

fH{length(fH)+1} = figure; hold on
    for whatMu=[1,3]
        for whatGass=1:length(gassType)+1 % +1 for LMS
            if whatMu==1
                whatLabel = strcat(gassTypeLabels{whatGass},' : ');
            else
                whatLabel = '';
            end
            plot( err.weight{whatMu,whatGass}(2,:) , ...    % error of weight 2 - i.e. a(2)
                'DisplayName', sprintf('%s $\\mu=%.2f$',whatLabel, mu(whatMu)), ...
                'Color', COLORS(whatGass,:), 'LineStyle', lineStyleList{whatMu}  );

        end
    end

    title('Weight Error with Time');
    xlabel('Time');
    ylabel('Weight Value', 'Color', 'k');
    grid minor;
    yLims = ylim;
    ylim([-inf, yLims(2)])
    xlim([0 600])
    legend('show','Orientation', 'vertical', 'NumColumns',2,'Location','best','FontSize',9)

    
lineStyleList = {'-','--',':'};
    
fH{length(fH)+1} = figure; hold on
    for whatMu=[1,3]
        for whatGass=1:length(gassType)+1 % +1 for LMS
            if whatMu==1
                whatLabel = strcat(gassTypeLabels{whatGass},' : ');
            else
                whatLabel = '';
            end
            plot( squeeze(pow2db( mean(err.diff{whatMu,whatGass}.^2, 1) )) , ...    % mean square error
                'DisplayName', sprintf('%s $\\mu=%.2f$',whatLabel, mu(whatMu)), ...
                'Color', COLORS(whatGass,:), 'LineStyle', lineStyleList{whatMu}  );

        end
    end

    title(' Mean Square Error with Time');
    xlabel('Time');
    ylabel('Error Value', 'Color', 'k');
    grid minor;
    yLims = ylim;
    ylim([-inf, yLims(2)*1.2])
    xlim([0 800])
    legend('show','Orientation', 'vertical', 'NumColumns',2,'Location','best','FontSize',9)
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q2_2a_fig',num2str(ii,'%02i')],'pdf')
    end
end
