function [X, y] = arima2diffEqns(signal, lags, varargin)
% arima2diffEqns	Convert time series to difference equation matrix.
% Input: 
%       - s: input signal,   [1 N] or [N 1]
%       - lags: time delays, [1 M] or [M 1] length(lags) = M (Order of Filter)
%               A lag of 0 is the current time-sample, for predictors you
%               want to start with a lag of 1, this could also be done by
%               specifying a Delta of 1
%       varargin:                                  variable input arguments
%       - varargin{1} = Delta: additional lag delay 
% Output: 
%       * X: Design matrix, rows is the signal delayed by lags [M N]
%       * y: Target vector with maximum delay,                 [1 N]
%   [X, y] = arima2diffEqns(s, lags, Delta) splits s into difference equations,
%            with lags specified by lags, and Delta. Note: the order of lags
%            in the lag vector matches the lag order of the output matrix

    Delta = 0;
    
    % check if input signal is 1D
    if ~isvector(signal)
        error("Input signal must be 1D");
    end
    % check if lags is vector of delays
    if ~isvector(lags)
        error("Lags parameter must be a 1D vector of delays, 0 Delay = Current Sample");
    end
    
    if ~isempty(varargin)
        if isa(uint64(varargin{1}),'integer')
            Delta = varargin{1};
            if Delta < 0
                error('Additional lag (Delta) should be positive')
            end
        else
            error('Additional lag (Delta) should be an integer')
        end
    end
    
    N = length(signal); % Length of the signal

    % Design matrix: pre-allocate for speed
    X = zeros(length(lags), N);
    % Target vector: pre-allocate for speed
    y = zeros(1, N);
    
    % update lags according to Delta
    lags = lags + Delta; 
    
    % Convolve signal with mask to get delay differences
    for nLag=1:length(lags)
        % Feature Matrix, each row is differently delayed
        X(nLag, :) = conv( [ zeros(lags(nLag), 1) ; 1 ], ...
                           signal(1:end-lags(nLag)) );
    end
    
    % Target Vector y
    y = conv( [ zeros(max(lags), 1) ; 1 ], ...
                signal(1:end-max(lags)) );
                          
end
