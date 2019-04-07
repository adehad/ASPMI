%% Q1 Classical and Modern Spectrum Estimation
% 1.6b
%{
Using only the r most significant principal components 
(as determined by the identified rank), create a low-rank approximation 
of Xnoise, denoted by X˜noise. 

Compare the difference (error) between the variables (columns) of the
noiseless input matrix, X, and those in the noise corrupted matrix Xnoise 
and denoised matrix X˜noise.
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

%% SVD Analysis

% SVD - matrix decomposition     U*S*(V') = A
[U.XNoise, S.XNoise, V.XNoise] = svd(Xnoise);

msErr = @(X1,X2) sum(sum(abs(X1 - X2).^2))/numel(X1); % MSE Lambda

rankComp = 1:6; % rank components to preserve


for ii=1:length(rankComp)
    Xnoise_hat =  U.XNoise(  :            , 1:rankComp(ii)) ...
                 *S.XNoise( 1:rankComp(ii), 1:rankComp(ii)) ...
                 *V.XNoise(  :            , 1:rankComp(ii)   )';
    
    estimateErr(ii) = msErr(X,Xnoise_hat); % MSE
%     estimateErr(ii) = norm(Xnoise_hat-X, 'fro');
end
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

fH{length(fH)+1} = figure;
    plot(rankComp, estimateErr);
    title('Low Rank Approximation Error');
    xlabel('Dimensionality Reduction Rank ');
    ylabel('Mean Squared Error');
    xticks(rankComp)
    grid minor;
    

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_6b_fig',num2str(ii,'%02i')],'pdf')
    end
end
