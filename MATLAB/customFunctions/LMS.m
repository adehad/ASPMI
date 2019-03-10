function [y, e, W] = LMS(X, d, mu, varargin)
% LMS	Least Mean Square (LMS) adaptive filter.
%       - X: design matrix, size(X)=[M N]
%       - d: target vector, size(d)=[1 N]
%       - mu: step size, scalar
%       - varargin(1) = gamma: leakage coefficient, scalar
%       * y: filter output, size(y)=[1 N]
%       * e: prediction error, d(n) - y(n)
%       * W: filter weights, size(W)=[M N]
%   [y, e, w] = LMS(X, d, mu, gamma) train LMS filter on Xd data.
    
    % check if design matrix is 2D
    if ~ismatrix(X)
        error("design matrix must be 2D");
    end
    % check if target vector is 1D
    if ~ismatrix(d)
        error("target vector must be 1D");
    end
    % check X-d size consistency
    if size(X, 2) ~= size(d, 2)
        error("design matrix and target vector sizes are inconsistent");
    end
    % check if step-size is scalar
    if ~isscalar(mu)
        error("step-size parameter must be scalar");
    end
    % check if leakage coefficient is scalar
    if ~isempty(varargin)
        gamma = varargin{1};
        if ~isscalar(gamma)
            error("leakage coefficient parameter must be scalar");
        end
        if (gamma>1 && gamma<0)
            warning('Gamma provided, %d, not between 0 and 1. \n Setting gamma to zero ', gamma)
            gamma = 0;
        end
    else
        warning('No leakage coefficient specified, using default: 0')
        gamma = 0;
    end
    
    % sizes
    [M, N] = size(X);
    % filter output: init
    y = zeros(size(d));
    % prediction error: init
    e = zeros(size(d));
    % LMS filter weights: init
    W = zeros(M, N);
    
    % iterate over time
    for n=1:N
        % filter output n, y(n)
        y(n) = W(:, n)' * X(:, n);
        % prediction error n, e(n)
        e(n) = d(n) - y(n);
        % weights update rule
        W(:, n+1) = (1 - mu * gamma) * W(:, n) + mu * e(n) * X(:, n);
    end
    
    % discard first weight
    W = W(:, 2:end);
end