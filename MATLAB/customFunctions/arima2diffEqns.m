function [X, y] = arima2diffEqns(signal, lags, varargin)
% XY	Convert time series to supervised learning (xi, yi) pairs, for ALE.
%       - s: input signal, size(s)=[1 N] or size(s)=[N 1]
%       - lags: time delays, 1D vector - length(lags) = M - Order of Filter
%       - varargin{1} = Delta: additional output delay 
%       * X: design matrix, size(X)=[M N-1]
%       * y: target vector, size(d)=[1 N-1]
%   [X, y] = XyDelta(s, Delta, M) splits s to AR(M) features-target pairs,
%   with Delta delays.

    % check if input signal is 1D
    if ~isvector(signal)
        error("input signal must be 1D");
    end
    % check if lags is vector of delays
    if ~isvector(lags)
        error("lags parameter must be a 1D vector of delays");
    end
    
    if ~isempty(varargin)
        Delta = varargin{1};
        if Delta < 0
            error('Additional Output lag (Delta) should be positive')
        end
    else
        Delta = 0;
    end
    
    % turn to column vector
%     signal = reshape(signal, [], 1);
    N = length(signal); % Length of the signal

    % design matrix: init
    X = zeros(length(lags), N);
    % target vector: init
    y = zeros(1, N);
    
    % update lags according to Delta
    lags = lags + Delta; 
    
    % Convolve signal with mask to get delay differences
    for nLag=1:length(lags)
        % feature vector n, x(n)
        X(nLag, :) = conv( [ zeros(lags(nLag), 1) ; 1 ], ...
                           signal(1:end-lags(nLag)) );
    end
    
    % target value n, y(n)
    y = conv( [ zeros(max(lags), 1) ; 1 ], ...
                signal(1:end-max(lags)) );
                          
end
