%% Q1 Classical and Modern Spectrum Estimation
%  1.2b
%{
Apply the standard periodogram approach to the entire recording,
as well as the averaged periodogram with different window lengths:
 (10 s, 5 s, 1 s) to the EEG data.
 Can you identify the the peaks in the spectrum corresponding to SSVEP? 
There should be a peak at:
 - the same frequency as the frequency of the ?ashing stimulus (integer X in
the range [11, . . . , 20]), known as the fundamental frequency response peak,
 - and at some integer multiples of this value, known as the harmonics of the response.
 It is important to note that the subject was tired during the recording
which induced a strong response within 8-10 Hz (so called alpha-rhythm), 
this is not the SSVEP. Also note that a power-line interference was induced 
in the recording apparatus at 50 Hz, and this too is not the SSVEP.
 
To enable a fair comparison across all spectral analysis approaches,
 you should keep the number of frequency bins the same.
**** Hint: It is recommended to have 5 DFT samples per Hz. **** 
%}
%% Premable
% Use Ctrl+Enter to run code section by section

% Clear out the workspace
clc; clear variables;
% 
% % Initialise common functions for this question
questionNum = 1;
q_initialise;

% Set the SAVE_FIGS to true if you want to save all figures
SAVE_FIGS = true;


%% Initialising Input Data and Variables

EEG = load(['resources',filesep,'EEG_Data',filesep,'EEG_Data_Assignment1.mat']);
% EEG.fs: Sampling rate in Hz = 1200 Hz
% EEG.POz: EEG in Volts, 80 seconds. 
    % incl. X*n Hz (stimulus + harmonics), 8-10 Hz(alpha-rythm), 50 Hz Mains noise
    
% Recenter lambda function (removes mean bias)
recenter = @(x) (x - mean(x));

deltaF = 1/5; % 5 DFT samples per Hz

N = size(EEG.POz,1); % Total Signal Length
K = EEG.fs/deltaF; % 6000 - Number of DFT points for 5 DFT points per Hz
    % Seen here:
    % https://dsp.stackexchange.com/questions/17579/number-of-dft-fft-points-required-for-a-specific-frequency-resolution-for-an-o
%% Initialising a Window
% Requires Signal Processing ToolBox
% https://uk.mathworks.com/help/signal/windows.html

% Run: windowDesigner
% for an interactive look at window types

win_BKH = blackmanharris(N); % Blackman Harris Window
winLabel = ["BK-Harris"];

% win_

%% Periodogram / PSD Calculations
legendString = [];

% Standard - Window = length of signal
[P_EEG.stan, fAx.stan, P_EEG_Conf.stan] = ... 
    periodogram(EEG.POz, rectwin(N), K, EEG.fs, 'onesided', 'power');
legendString = [legendString, "Standard Window"];

% Averaged Periodogram
% Note: We use the Welch estimate and specify an overlap of 0 to achieve
%       the Bartlett estimate
overlap = 0; % samples of window overlap, if not specified defaults to 50% in pwelch
% Window = 10s
[P_EEG.win10s, fAx.win10s, P_EEG_Conf.win10s] = ...
    pwelch(EEG.POz, rectwin(10*EEG.fs), overlap, K, EEG.fs, 'onesided', 'power','ConfidenceLevel',0.95);
legendString = [legendString, "10s Window"];

% Window = 5s
[P_EEG.win5s, fAx.win5s, P_EEG_Conf.win5s] = ...
    pwelch(EEG.POz, rectwin(5*EEG.fs), overlap, K, EEG.fs, 'onesided', 'power','ConfidenceLevel',0.95);
legendString = [legendString, "{ }5s Window"];

% Window = 1s
[P_EEG.win1s, fAx.win1s, P_EEG_Conf.win1s] = ...
    pwelch(EEG.POz, rectwin(1*EEG.fs), overlap, K, EEG.fs, 'onesided', 'power','ConfidenceLevel',0.95);
legendString = [legendString, "{ }1s Window"];

% NOTE: All the fAx are the same as we have fixed the number of DFT
%       points to 5 per Hz

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
% fAx = 0:EEG.fs/K:1-EEG.fs/K;

fH{1} = figure;
    plot(fAx.stan,pow2db(P_EEG.stan))
    title("Full Size Rectangular Window EEG Periodogram");
    xlabel("Frequency (Hz)");
    ylabel("Power Density (dB)");
    grid minor;
    xlim([0 60])
    
fH{2} = figure;
    hold on
    plotH(1) = plot(fAx.win10s,pow2db(P_EEG.win10s));
    plotH(2) = plot(fAx.win5s,pow2db(P_EEG.win5s));
    plotH(3) = plot(fAx.win1s,pow2db(P_EEG.win1s));
    
    % Confidence Levels
    plot(fAx.win10s,pow2db(P_EEG_Conf.win10s),'Color',COLORS(1,:),'LineWidth',0.5,'LineStyle','-.')
    plot(fAx.win5s,pow2db(P_EEG_Conf.win5s),'Color',COLORS(2,:),'LineWidth',0.5,'LineStyle','-.')
    plot(fAx.win1s,pow2db(P_EEG_Conf.win1s),'Color',COLORS(3,:),'LineWidth',0.5,'LineStyle','-.')

    %
        
    title("Varied Window Size EEG Bartlett Average Periodogram");
    xlabel("Frequency (Hz)");
    ylabel("Power Density (dB)");
    grid minor;
    xlim([0 60])
    legend( plotH, legendString([2:end]), 'NumColumns', 3, 'Location', 'South' )

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q1_2b_fig',num2str(ii,'%02i')],'pdf')
    end
end
