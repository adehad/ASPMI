%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
% 3.1c
%{
Generate two sets of complex voltages, one balanced and one unbalanced.
 To generate an unbalanced system change the magnitude and/or phase of one
 or more phases. 
Plot the circularity diagrams of these complex alpha - beta voltages. 

Comment on the shape of the circularity diagram when the system is
balanced vs. unbalanced.

How would you use the circularity diagram to identify a fault in the system?
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
N = 1e3;
% mains frequency
f0 = 50;
% sampling rate
fs = 50e3; 

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

%% Realisations 

clarkeVoltage.balanced   = clarkeVec(clarkeTrans( voltageGen(V.balanced,Delta.balanced) ));
clarkeVoltage.unbalanced = clarkeVec(clarkeTrans( voltageGen(V.unbalanced,Delta.unbalanced) ));

dDelta=0.25:0.25:1;
for ii=1:length(dDelta) % keep voltage balanced and sweep delta
    newDelta = Delta.balanced + [ 0; -dDelta(ii); +dDelta(ii) ];
    clarkeVoltage.dDel(ii,:) = clarkeVec(clarkeTrans( voltageGen(V.balanced,newDelta) ));
end

dVoltage=0.25:0.25:1;
for ii=1:length(dVoltage) % keep delta balanced and sweep voltage
    newVoltage = V.balanced + [ 0; -dVoltage(ii); +dVoltage(ii) ];
    clarkeVoltage.dVolt(ii,:) = clarkeVec(clarkeTrans( voltageGen(newVoltage,Delta.balanced) ));
end

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

balanceTypes = fieldnames(clarkeVoltage);
balanceLabels = {"Balanced","Unbalanced"};

% sample circularity investigation
for ii=1:2
    tempClarke = clarkeVoltage.(balanceTypes{ii});
    fH{length(fH)+1} = figure; hold on
        scatter(real(tempClarke), imag(tempClarke), 30, COLORS(ii,:), 'filled')
        title( sprintf('%s, $|\\rho|$=%.3f \n V=[%.2f;%.2f;%.2f], $\\Delta$=[%.2f,%.2f,%.2f]',...
                balanceLabels{ii}, circularity(tempClarke), V.(balanceTypes{ii})(1), V.(balanceTypes{ii})(2), V.(balanceTypes{ii})(3),...
                                                            Delta.(balanceTypes{ii})(1), Delta.(balanceTypes{ii})(2), Delta.(balanceTypes{ii})(3)))
        xlabel('Real Part, $\Re$');
        ylabel('Imaginary Part, $\Im$');

        grid minor
end



% phase & amplitude sweeps
sweepType = {' ',' ','Phase','Magnitude'};
sweepSymb = {' ',' ','\Delta','V'};
for ii=3:4
    fH{length(fH)+1} = figure; hold on
    
        for jj=1:size(clarkeVoltage.(balanceTypes{ii}),1)
            tempClarke = clarkeVoltage.(balanceTypes{ii})(jj,:);
            if strcmpi(balanceTypes{ii},'dDel'); plotLabels = dDelta; else; plotLabels = dVoltage; end

            scatter(real(tempClarke), imag(tempClarke), ...
                    20, COLORS(jj,:), 'filled', ...
                    'DisplayName', sprintf('$%s$: %.2f',sweepSymb{ii},plotLabels(jj)))
        end
        title(sprintf('%s Parameter Sweep',sweepType{ii}))
        xLims = xlim;
%         xlim([xLims(1), xLims(2)*1.5])
        xlabel('Real Part, $\Re$');
        ylabel('Imaginary Part, $\Im$');
        legend('show','Location','best','NumColumns',length(plotLabels))
        grid minor
        ylim([-1.8,2])
        xlim([-1.75,1.75])

end 

%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_1c_fig',num2str(ii,'%02i')],'pdf')
    end
end
