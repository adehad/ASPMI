function [y, e, W] = NLMS(X, d, mu, varargin)
% NLMS	Normalised Least Mean Square (NLMS) adaptive filter.
% Input: 
%       - X: Design matrix, [M N]
%       - d: Target vector, [1 N]
%       - mu: Step size, numeric
%       varargin:                                  variable input arguments
%       - varargin{1}   = setGNGD: GNGD Flag (default: false)
%       - varargin{2}   = rho:               (default: 0)
% Output: 
%       * y: Filter output,     [1 N]
%       * e: Prediction error,  d-y
%       * W: Filter weights,    [M N]
% Usage: 
%   [y, e, w] = LMS(X, d, mu, setGNGD, rho) train NLMS filter on Xd data.
    
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
    
    % Check if GNGD flag set, if it is, get the rho value
    setGNGD = false;
    if ~isempty(varargin)
        if strcmpi(varargin{1},'gngd')
            setGNGD = true;
            if length(varargin)>1
                rho = varargin{2};
            else
                error('must specify rho value when using GNGD')
            end
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
    % Regularization Factor: 1 / mu
%     epsilon = ones(N+1, 1) * eps; % / mu;
    epsilon = ones(N+1, 1) / mu;
    
    beta = 1;
    
    % Iterate over the discrete time samples
    for n=1:N
        % Filter output n, y(n)
        y(n) = W(:,n)' * X(:,n);
        % Prediction error n, e(n)
        e(n) = d(n) - y(n);
        % Weights update rule
        W(:,n+1) = W(:,n) + ...
                   beta*e(n)*X(:,n)  ...
                       /( epsilon(n) + X(:,n)'*X(:,n) ) ;
        if n > 1
            if setGNGD
            % epsilon update rule
            epsilon(n+1) = epsilon(n) - ...
                            rho*mu* ( e(n)*e(n-1) * X(:, n)'*X(:, n-1) )  ...
                                  / ( epsilon(n-1) + X(:, n-1)'*X(:, n-1) )^2;
            end
        end
    end
    
    % Discard first weight
    W = W(:,2:end);
    
    % Check Instability
    if find(isnan(y)==1,1)
        warning('unstable mu provided, output reached NaN')
    end
end