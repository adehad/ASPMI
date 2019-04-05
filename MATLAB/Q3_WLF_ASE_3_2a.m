%% Q3 Widely Linear Filtering and Adaptive Spectrum Estimation
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
N = 1.5e3;
% sampling rate
fs = 1.5e3; 

% noise
sigma_sq = 0.05;
eta = wgn(N, 1, pow2db(sigma_sq), 'complex');


% make frequencies
freq = zeros(1,N);  % pre-allocate
thresh = [5e2,1e3,1.5e3];        % thresholds
for n=1:length(freq)
    if n <= thresh(1)
        freq(n) = 100; 
    elseif n <= thresh(2)
        freq(n) = 100 + (n-500)/2;
    elseif n <= thresh(3)
        freq(n) = 100 + ( (n-1000)/25 )^2;
    else
        freq(n) = 0;
    end
end

% make phase - by integration
phi = cumtrapz(freq);

msErr = @(x) mean(abs(x).^2);
%% AR 
% noisy signal
y = exp( (2*pi* phi/fs)* 1j ) + eta.';

ar_model{1}.p = [1 5 10];
for ii=1:length(ar_model{1}.p)
    % AR(p)
    a = aryule( y, ar_model{1}.p(ii) );
    
    % power density estimate
    [ar_model{1}.h(ii,:), ar_model{1}.w(ii,:)] = freqz(1, a, N, fs);
    
    % psd
    ar_model{1}.psd(ii,:) = mag2db( abs(ar_model{1}.h(ii,:)) );
end


ar_model{2}.p = 1;
ar_model{2}.thresh = [1, thresh];
for jj=1:length(ar_model{2}.thresh)-1
    for ii=1:length(ar_model{2}.p)
        % segment to look at
        seg = ar_model{2}.thresh(jj):ar_model{2}.thresh(jj+1);
        
        % AR(p)
        a = aryule( y(seg), ar_model{2}.p );

        % power density estimate
        [ar_model{2}.h(jj,:,ii), ar_model{2}.w(jj,:,ii)] = freqz(1, a, N, fs);

        % psd
        ar_model{2}.psd(jj,:,ii) = mag2db( abs(ar_model{2}.h(jj,:,ii)) );
    end
end

%% Plots
close all % close current figures
fH = []; % clear the figure handle variable
plotH = []; % clear the plot handle variable
legendString = []; % clear the legend string variable

% plot frequencies
fH{length(fH)+1} = figure; hold on
%     yyaxis left
    plot(1:N, freq)
    title("FM: Frequency, $f(n)$");
    xlabel("Time Index, $n$");
    ylabel("Frequency (Hz)");
    grid on; grid minor;
    yLim = ylim;
    ylim([0,yLim(2)])
    
% plot phase
fH{length(fH)+1} = figure; hold on
%     yyaxis right
    subplot(2,1,1)
        plot( 1:N, rad2deg(wrapToPi(phi)), 'Color', COLORS(2,:) )
        title("FM: Phase, $\phi(n)$");
    %     xlabel("Time Index, $n$");
        ylabel("Phase ($^{\circ}$)");
        grid on; grid minor;
        xlim([250, 750])
        yticks([-360:90:360])
        set(gca,'Position',POSITION.subplot211)
    subplot(2,1,2)
        plot( 1:N, rad2deg(wrapToPi(phi)), 'Color', COLORS(2,:) )
        xlabel("Time Index, $n$");
        ylabel("Phase ($^{\circ}$)");
        grid on; grid minor;
        xlim([750, 1250])
        yticks([-360:90:360])
        set(gca,'Position',POSITION.subplot212)
     
        
fH{length(fH)+1} = figure; hold on
    for ii=1:length(ar_model{1}.p)
        plot(ar_model{1}.w(ii,:), ar_model{1}.psd(ii,:),'DisplayName', sprintf('$p$=%i',ar_model{1}.p(ii)))
    end
    title("FM: \texttt{aryule}-AR");
    xlabel("Frequency (Hz)");
    ylabel("PSD (db)");
    grid on; grid minor;
    legend('show')
    xlim([0, 700])
    
fH{length(fH)+1} = figure; hold on
    for jj=1:length(ar_model{2}.thresh)-1
        for ii=1:length(ar_model{2}.p)
            plot(ar_model{2}.w(jj,:,ii), ar_model{2}.psd(jj,:,ii),'DisplayName', sprintf('%i:%i',ar_model{2}.thresh(jj),ar_model{2}.thresh(jj+1)))
        end
    end
    title(sprintf("FM: \texttt{aryule}-AR(%i) ",ar_model{2}.p(ii)));
    xlabel("Frequency (Hz)");
    ylabel("PSD (db)");
    grid on; grid minor;
    legend('show')
    xlim([0, 700])
    
%% Save Figures

if SAVE_FIGS
    for ii=1:length(fH) % For all figure handles
        saveas(fH{ii},['figures', filesep,'q3_2a_fig',num2str(ii,'%02i')],'pdf')
    end
end
