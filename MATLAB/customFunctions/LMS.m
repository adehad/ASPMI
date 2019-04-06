function [y, e, W] = LMS(X, d, mu, varargin)
% LMS	Least Mean Square (LMS) adaptive filter.
% Input: 
%       - X: Design matrix, [M N]
%       - d: Target vector, [1 N]
%       - mu: Step size, numeric
%       varargin:                                  variable input arguments
%       - varargin{1}   = gamma: leakage coefficient, scalar (default: 0)
% Output: 
%       * y: Filter output,     [1 N]
%       * e: Prediction error,  d-y
%       * W: Filter weights,    [M N]
% Usage: 
%   [y, e, w] = LMS(X, d, mu, gamma) train LMS filter on Xd data.
    
    gamma = 0;

    % Design matrix is 2D
    if ~ismatrix(X)
        error("Design matrix must be 2D, [M N]");
    end
    
    % Target / Ground Truth is 1D
    if ~isvector(d)
        error("Target vector must be 1D, [1 N]");
    end
    
    % X-d Size Match
    if size(X, 2) ~= size(d, 2)
        if size(X, 2) == size(d.', 2)
            d = d.';    % Using MATLAB {.'} operator to prevent conjugate transpose of complex data
            warning('Auto-transposing target matrix data')
        else
            error("Design matrix and target vector sizes are incompatible, [M N] and [1 N] required");
        end
    end
    
    % Step-size is a numeric scalar
    if ~isa(mu,'numeric')
        error("Step-size parameter (mu) must be numeric");
    end
    
    % Check if leakage coefficient is numeric
    if ~isempty(varargin)
        if isa(varargin{1},'numeric')
            gamma = varargin{1};
            if (gamma>1 && gamma<0)
                warning('Gamma provided, %d, not between 0 and 1. \n Setting gamma to zero ', gamma)
            end
        else
            error("leakage coefficient parameter must be scalar");
        end
    end
    
    % sizes
    [M, N] = size(X);
    % Filter Output: pre-allocate for speed
    y = zeros(size(d));
    % Prediction Error: pre-allocate for speed
    e = zeros(size(d));
    % LMS filter weights: pre-allocate for speed
    W = zeros(M, N);
    
    % Iterate over the discrete time samples
    for n=1:N
        % Filter output n, y(n)
        y(n) = W(:,n)' * X(:,n);
        % Prediction error n, e(n)
        e(n) = d(n) - y(n);
        % Weights update rule
        W(:,n+1) = (1 - mu*gamma) * W(:,n) + mu * e(n) * X(:,n);
    end
    
    % Discard first weight
    W = W(:,2:end);
    
    % Check Instability
    if find(isnan(y)==1,1)
        warning('unstable mu provided, output reached NaN')
    end
end