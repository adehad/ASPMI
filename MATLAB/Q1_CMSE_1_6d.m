%% Q1 Classical and Modern Spectrum Estimation
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;

% Initialise common functions for this question
questionNum = 1;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
SAVE_FIGS = true;


%% LOAD PCA/PCR DATA
load(['resources',filesep,'PCR',filesep,'PCAPCR.mat']);
addpath(['resources',filesep,'PCR',filesep]); % contains regval()
% X
% Xnoise
% Xtest
% Y
% Ytest

%% PCR and OLS Analysis
% msErr = @(X1,X2) sum(sum(abs(X1 - X2).^2))/numel(X1);
msErr = @(X1,X2) norm(abs(X1 - X2),2).^2; % 2-norm squared

% PCR 
% SVD - matrix decomposition     U*S*(V') = A
[U.XNoise, S.XNoise, V.XNoise] = svd(Xnoise);

rankComp = 3:8; % rank components to preserve

% OLS
% weights
OLS.B = Xnoise' * Xnoise \ Xnoise' * Y; % \ performs the inverse


% Generate Realisations (OLS)
numRealisations = 100;
for ii=1:numRealisations
    % Simulation
    [OLS.Y_hat, OLS.Y] = regval(OLS.B);
    % Store Error
    ERRS.OLS(ii) = msErr(OLS.Y_hat,OLS.Y);
end

for ii=1:length(rankComp)
    
    % X = U*S*(V')
    Xnoise_hat =  U.XNoise(  :            , 1:rankComp(ii)) ...
                 *S.XNoise( 1:rankComp(ii), 1:rankComp(ii)) ...
                 *V.XNoise(  :            , 1:rankComp(ii)   )';
     
             
    % weights - B_PCR = V*S*(U')*Y
    PCR.B =   V.XNoise(  :            , 1:rankComp(ii)) ...
             /S.XNoise( 1:rankComp(ii), 1:rankComp(ii)) ... % \ performs the inverse
             *U.XNoise(  :            , 1:rankComp(ii))' ...
             *Y;
         
    for jj=1:numRealisations
        % Simulation
        [PCR.Y_hat, PCR.Y] = regval(PCR.B);
        % Store Error
        ERRS.PCR(ii,jj) = msErr(PCR.Y_hat, PCR.Y);
    end
end
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

fH{length(fH)+1} = figure;
    plot([rankComp(1) rankComp(end)], [ERRS.OLS' ERRS.OLS'], 'Color', COLORS(6, :),'LineWidth',0.5);
    hold on
    plot(rankComp, ERRS.PCR, 'Color', COLORS(3, :),'LineWidth',0.5);
    % Mean Trace
    plotH(1)= plot([rankComp(1) rankComp(end)], [mean(ERRS.OLS)' mean(ERRS.OLS)'], 'Color', COLORS(1, :),'DisplayName',"OLS");
    plotH(2)= plot(rankComp, mean(ERRS.PCR,2), 'Color', COLORS(2, :),'DisplayName',"PCR");
    hold off
    
    title("OLS vs PCR Estimation Error");
    xlabel("Dimensionality Reduction Rank ");
    ylabel("Mean Squared Error");
    xticks(rankComp)
    grid minor;
    legend(plotH)

    

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_6d_fig',num2str(ii,'%02i')],'pdf')
    end
end
