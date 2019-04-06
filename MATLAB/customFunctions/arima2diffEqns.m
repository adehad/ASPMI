function [X, y] = arima2diffEqns(signal, lags, varargin)
% arima2diffEqns	Convert time series to differential equation matrix.
% Input: 
%       - s: input signal,  [1 N] or [N 1]
%       - lags: time delays, 1D vector - length(lags) = M (Order of Filter)
%       varargin:                                  variable input arguments
%       - varargin{1} = Delta: additional output delay 
% Output: 
%       * X: Design matrix, [M N-1]
%       * y: Target vector, [1 N-1]
%   [X, y] = arima2diffEqns(s, lags, Delta) splits s differential equations,
%            with lags specified by lags, and Delta.

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
        if isa(varargin{1},'integer')
            Delta = varargin{1};
            if Delta < 0
                error('Additional Output lag (Delta) should be positive')
            end
        else
            error('Additional Output lag (Delta) should be an integer')
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
