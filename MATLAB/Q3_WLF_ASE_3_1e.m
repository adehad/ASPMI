%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
% 3.1e
%{
Use the CLMS given (35) and ACLMS algorithms given in (36) to estimate 
the frequency of the alpha - beta voltages you generated in Part b).
 For unbalanced system voltages, does the CLMS give the correct frequency estimate?
 If not, why?
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

% number of samples
N = 2e3;
% mains frequency
f0 = 50;
% sampling rate
fs = 5e3; 

% voltage phase
phi = [ 0;
        -2*pi/3;
        +2*pi/3 ];

% additional phase shift
Delta.balanced   = zeros(size(phi));
Delta.unbalanced = [0.5; 3; 5];         % we change phase and/or magnitude

% amplitude
V.balanced   = ones(size(phi));
V.unbalanced = randn(size(phi))*3;      % we change phase and/or magnitude


voltageGen = @(V,Delta) V.*cos(2*pi*( f0/fs )*(1:N) +Delta +phi);

% Lecture 6, Slide 11:
circularity = @(data) abs( mean(data*data.') /mean(data*data') ); % NOTE: {.'} vs {'} operators
% I would have thought a better one would be to find major and minor axis
% by PCA and compare magnitude

% Define the Clarke Transform, aka: alpha-beta transform
clarkeTrans = @(voltage)  sqrt(2/3)*[ sqrt(2)/2, sqrt(2)/2  ,    sqrt(2)/2  ; ...
                                         1     ,    -1/2    ,       -1/2    ; ... Clarke Matrix
                                         0     , sqrt(3)/2  ,    -sqrt(3)/2 ] ... 
                                   *voltage; % Apply to input voltage vector
                      
clarkeVec = @(clarkeT) complex( clarkeT(2,:), clarkeT(3,:) );  

msErr = @(x) mean(abs(x).^2);

% step-size 
mu =  0.04; 
% model orders
M = 1;
lags = @(order) 0:order-1; % starts at 0 lag, hence -1
%% Realisations 

clarkeVoltage.balanced   = clarkeVec(clarkeTrans( voltageGen(V.balanced,Delta.balanced) ));
clarkeVoltage.unbalanced = clarkeVec(clarkeTrans( voltageGen(V.unbalanced,Delta.unbalanced) ));

% dDelta=0.25:0.25:1;
% for ii=1:length(dDelta) % keep voltage balanced and sweep delta
%     newDelta = Delta.balanced + [ 0; -dDelta(ii); +dDelta(ii) ];
%     clarkeVoltage.dDel(ii,:) = clarkeVec(clarkeTrans( voltageGen(V.balanced,newDelta) ));
% end
% 
% dVoltage=0.25:0.25:1;
% for ii=1:length(dVoltage) % keep delta balanced and sweep voltage
%     newVoltage = V.balanced + [ 0; -dVoltage(ii); +dVoltage(ii) ];
%     clarkeVoltage.dVolt(ii,:) = clarkeVec(clarkeTrans( voltageGen(newVoltage,Delta.balanced) ));
% end

%% CLMS
    %% Balanced
    % create differential equations
    [X, ~] = arima2diffEqns(clarkeVoltage.balanced, lags(M),1);

    % CLMS
    [~, err.CLMS{1}, weight.H_CLMS{1}] = CLMS(X, clarkeVoltage.balanced, mu);
    % ACLMS
    [~, err.ACLMS{1}, weight.H_ACLMS{1}, weight.G_ACLMS{1}] = CLMS(X, clarkeVoltage.balanced, mu, 'aug');

    % f0 estimates
    est_f0.CLMS{1}  = fs/(2*pi) * atan( imag(weight.H_CLMS{1}) ./ real(weight.H_CLMS{1}) );
    est_f0.ACLMS{1} = fs/(2*pi) * atan( sqrt( imag(weight.H_ACLMS{1}).^2 - abs(weight.G_ACLMS{1}).^2 ) ./ real(weight.H_ACLMS{1}) );
    %% UnBalanced
    % create differential equations
    [X, ~] = arima2diffEqns(clarkeVoltage.unbalanced, lags(M),1);

    % CLMS
    [~, err.CLMS{2}, weight.H_CLMS{2}] = CLMS(X, clarkeVoltage.unbalanced, mu);
    % ACLMS
    [~, err.ACLMS{2}, weight.H_ACLMS{2}, weight.G_ACLMS{2}] = CLMS(X, clarkeVoltage.unbalanced, mu, 'aug');

    % f0 estimates
    est_f0.CLMS{2}  = fs/(2*pi) * atan( imag(weight.H_CLMS{2}) ./ real(weight.H_CLMS{2}) );
    est_f0.ACLMS{2} = fs/(2*pi) * atan( sqrt( imag(weight.H_ACLMS{2}).^2 - abs(weight.G_ACLMS{2}).^2 ) ./ real(weight.H_ACLMS{2}) );
%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

balanceTypes = fieldnames(clarkeVoltage);
balanceLabels = {"Balanced","Unbalanced"};

t_steady = 800;

for ii=1:2
    %% error plot
    tempClarke = clarkeVoltage.(balanceTypes{ii});
    fH{length(fH)+1} = figure; hold on
        plot(pow2db(abs(err.CLMS{ii}).^2), 'DisplayName', 'CLMS')
        plot(pow2db(abs(err.ACLMS{ii}).^2), 'DisplayName', 'ACLMS')
        title( sprintf('%s Error \n V=[%.2f;%.2f;%.2f], $\\Delta$=[%.2f,%.2f,%.2f]',...
                balanceLabels{ii},  V.(balanceTypes{ii})(1), V.(balanceTypes{ii})(2), V.(balanceTypes{ii})(3),...
                                    phi(1) + Delta.(balanceTypes{ii})(1), phi(2) + Delta.(balanceTypes{ii})(2), phi(3) + Delta.(balanceTypes{ii})(3)))
        xlabel('Time Index');
        ylabel('MSE (dB)');
        legend('show','Location','best')
        grid minor
    %% Frequency estimation 
    fH{length(fH)+1} = figure; hold on
        plot(abs(est_f0.CLMS{ii}), 'DisplayName', 'CLMS')
        plot(abs(est_f0.ACLMS{ii}), 'DisplayName', 'ACLMS')
        plot([0 N], [f0 f0], 'DisplayName', '$50\ Hz$', 'LineStyle', ':', 'Color', COLORS(6,:));
        title( sprintf('%s Frequency Estimate \n V=[%.2f;%.2f;%.2f], $\\Delta$=[%.2f,%.2f,%.2f]',...
                balanceLabels{ii},  V.(balanceTypes{ii})(1), V.(balanceTypes{ii})(2), V.(balanceTypes{ii})(3),...
                                    phi(1) + Delta.(balanceTypes{ii})(1), phi(2) + Delta.(balanceTypes{ii})(2), phi(3) + Delta.(balanceTypes{ii})(3)))
        xlabel('Time Index');
        ylabel('Frequency (Hz)');
        ylim([0 150])
        legend('show','Location','southeast')
        grid minor
end

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_1e_fig',num2str(ii,'%02i')],'pdf')
    end
end
