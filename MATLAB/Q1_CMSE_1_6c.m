%% Q1 Classical and Modern Spectrum Estimation
% 1.6c
%{
Calculate the OLS and PCR solutions for the parameter matrix B, 
which relates Xnoise and Y. 
Next, compare the estimation error between Y and YˆOLS = XnoiseBˆOLS and
                                                 YˆPCR = X˜noiseBˆPCR. 

Explain what happens when you estimate the data from the test-set using the
 regression coefficients computed from the training set, and quantify the
performance by comparing Ytest and Yˆtest-OLS = XtestBˆOLS with
                                   Yˆtest-PCR = X˜testBˆPCR.
%}
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
% X
% Xnoise
% Xtest
% Y
% Ytest

%% PCR and OLS Analysis
msErr = @(X1,X2) sum(sum(abs(X1 - X2).^2))/numel(X1);
% msErr = @(X1,X2) norm(abs(X1 - X2),2).^2; % 2-norm squared

% PCR 
% SVD - matrix decomposition     U*S*(V') = A
[U.XNoise, S.XNoise, V.XNoise] = svd(Xnoise);
[U.XTest, S.XTest, V.XTest] = svd(Xtest);

rankComp = 3:8; % rank components to preserve

% OLS
% weigts
OLS.B = Xnoise' * Xnoise \ Xnoise' * Y; % \ performs the inverse if a matrix
% targets
OLS.Y = Xnoise * OLS.B;
OLS.Y_test = Xtest * OLS.B;
% errors
ERRS.OLS_train = msErr(Y, OLS.Y);
ERRS.OLS_test  = msErr(Ytest, OLS.Y_test);
% ERRS.OLS_train = norm( Y - OLS.Y, 'fro');
% ERRS.OLS_test  = norm(Ytest - OLS.Y_test, 'fro');

for ii=1:length(rankComp)
    
    % X = U*S*(V')
    Xnoise_hat =  U.XNoise(  :            , 1:rankComp(ii)) ...
                 *S.XNoise( 1:rankComp(ii), 1:rankComp(ii)) ...
                 *V.XNoise(  :            , 1:rankComp(ii)   )';
    
    % X = U*S*(V')         
    Xtest_hat  =  U.XTest(  :            , 1:rankComp(ii)) ...
                 *S.XTest( 1:rankComp(ii), 1:rankComp(ii)) ...
                 *V.XTest(  :            , 1:rankComp(ii)   )';
    
             
    % weights - B_PCR = V*S^-1*(U')*Y
    PCR.B =   V.XNoise(  :            , 1:rankComp(ii)) ...
             /S.XNoise( 1:rankComp(ii), 1:rankComp(ii)) ... % / performs the inverse if a matrix
             *U.XNoise(  :            , 1:rankComp(ii))' ...
             *Y;
         
    % targets
    PCR.Y      = Xnoise_hat * PCR.B;
    PCR.Y_test = Xtest_hat  * PCR.B;
    % errors
    ERRS.PCR_train(ii) = msErr(Y, PCR.Y);
    ERRS.PCR_test(ii)  = msErr(Ytest, PCR.Y_test);
end
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

fH{length(fH)+1} = figure; hold on
    plot([rankComp(1) rankComp(end)], [ERRS.OLS_train ERRS.OLS_train],'DisplayName','OLS');
    plot(rankComp, ERRS.PCR_train,'DisplayName','PCR');
    hold off
    title('OLS vs PCR Error (Training)');
    xlabel('Dimensionality Reduction Rank ');
    ylabel('Mean Squared Error');
    xticks(rankComp)
    grid minor;
    legend('show')
    
fH{length(fH)+1} = figure; hold on
    plot([rankComp(1) rankComp(end)], [ERRS.OLS_test ERRS.OLS_test],'DisplayName','OLS');
    plot(rankComp, ERRS.PCR_test,'DisplayName','PCR');
    hold off
    title('OLS vs PCR Error (Testing)');
    xlabel('Dimensionality Reduction Rank ');
    ylabel('Mean Squared Error');
    xticks(rankComp)
    grid minor;
    legend('show')
    

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_6c_fig',num2str(ii,'%02i')],'pdf')
    end
end
