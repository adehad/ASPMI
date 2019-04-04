function [y, e, W] = NLMS(X, d, mu, varargin)
% NLMS	Normalised Least Mean Square (NLMS) adaptive filter.
%       - X: design matrix, size(X)=[M N]
%       - d: target vector, size(d)=[1 N]
%       - mu: step size, scalar
%       - varargin(1) = gamma: leakage coefficient, scalar
%       * y: filter output, size(y)=[1 N]
%       * e: prediction error, d(n) - y(n)
%       * W: filter weights, size(W)=[M N]
%   [y, e, w] = NLMS(X, d, mu, gamma) train LMS filter on Xd data.
    
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
    setGNGD = false;
    if ~isempty(varargin)
        if strcmpi(varargin{1},'gngd')
            setGNGD = true;
            rho = varargin{2};
        end
    end
    
    % sizes
    [M, N] = size(X);
    % filter output: init
    y = zeros(size(d));
    % prediction error: init
    e = zeros(size(d));
    % LMS filter weights: init
    W = zeros(M, N+1);
    % regularization factor: 1 / mu
%     epsilon = ones(N+1, 1) * eps; % / mu;
    epsilon = ones(N+1, 1) / mu;
    
    beta = 1;
    
    % iterate over time
    for n=1:N
        % filter output n, y(n)
        y(n) = W(:, n)' * X(:, n);
        % prediction error n, e(n)
        e(n) = d(n) - y(n);
        % weights update rule
        W(:, n+1) = W(:, n) + ...
                    beta*e(n)*X(:, n)  ...
                        /( epsilon(n) + X(:, n)'*X(:, n) ) ;
        if n > 1
            if setGNGD
            % epsilon update rule
            epsilon(n+1) = epsilon(n) - ...
                            rho*mu* ( e(n)*e(n-1) * X(:, n)'*X(:, n-1) )  ...
                                  / ( epsilon(n-1) + X(:, n-1)'*X(:, n-1) )^2;
            end
        end
    end
    
    % discard first weight
    W = W(:, 2:end);
end