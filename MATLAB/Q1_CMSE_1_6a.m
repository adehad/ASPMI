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
% X
% Xnoise
% Xtest
% Y
% Ytest

%% SVD Analysis

% SVD: Noise-free matrix
SVD.X = svd(X);
% SVD: Noise matrix
SVD.XNoise = svd(Xnoise);
% Squared error of singular values
sqError = (SVD.X - SVD.XNoise).^2;


%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

fH{length(fH)+1} = figure;
    bar([SVD.X,SVD.XNoise]);
    title("Singular Values of $X$ and $X_{noise}$");
    xlabel("Index");
    ylabel("SVD Magnitude");
    grid minor;
    legend(["pure","noisy"],'NumColumns',2)
    
fH{length(fH)+1} = figure;
    bar([sqError]);
    title("Square Error of $X$ and $X_{noise}$ Singular Values");
    xlabel("Index");
    ylabel("Square Error Magnitude");
    grid minor;

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_6a_fig',num2str(ii,'%02i')],'pdf')
    end
end
