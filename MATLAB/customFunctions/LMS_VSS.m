function [y, e, W] = LMS_VSS(X, d, mu_0, gamma, rho, gassType, varargin)
% GASS	Gradient Adaptive Step-Size (GASS) Least Mean Square (LMS) adaptive filter.
%       - X: design matrix, size(X)=[M N]
%       - d: target vector, size(d)=[1 N]
%       - mu_0: initial step size, scalar
%       - rho: learning rate, scalar
%       - gamma: leakage coefficient, scalar
%       - algo: algorithm name, string from
%               {'benveniste', 'ang_farhang', 'matthews_xie'}
%       - alpha: Ang & Farhang learning parameter, scalar
%       * y: filter output, size(y)=[1 N]
%       * e: prediction error, d(n) - y(n)
%       * W: filter weights, size(W)=[M N]
%   [y, e, W] = LMS_VSS(X, d, mu_0, rho, gamma, algo, alpha) train GASS filter on Xd data.

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
    % check if initial step-size is scalar
    if ~isscalar(mu_0)
        error("initial step-size parameter must be scalar");
    end
    % check if learning rate is scalar
    if ~isscalar(rho)
        error("learning rate parameter must be scalar");
    end
    % check if leakage coefficient is scalar
    if ~isscalar(gamma)
        error("leakage coefficient parameter must be scalar");
    end
    % check if algorithm type is string
    if ~isstring(gassType)
        error("algorithm type parameter must be string. e.g. ''benveniste'' OR ""benveniste"" ");
    end
    % check if Ang & Farhang learning parameter is set and is scalar
    if ~isempty(varargin)
        alpha = varargin{1};
        if ~isscalar(alpha)
            error("Ang & Farhang learning parameter (alpha) must be scalar");
        end
    else
        if strcmpi(gassType,"ang_farhang")
            error("Ang & Farhang learning parameter (alpha) must be specified");
        else
            alpha = 0; % if not using the Ang & Farhang method, we don't mind if it is zero
        end
    end

    % sizes
    [M, N] = size(X);
    % filter output: init
    y = zeros(size(d));
    % prediction error: init
    e = zeros(size(d));
    % GASS filter weights: init
    W = zeros(M, N+1);
    % step-size: init
    mu = zeros(1, N+1);
    mu(1) = mu_0;
    % phi term init:
    phi = zeros(M, N+1);
    
    % Establish what implementation of VSS to use
    if strcmpi(gassType, "benveniste")
        phiFunc = @(X,e,phi,mu,alpha) (eye(M) - mu*X*X')*phi + e*X;
        
    elseif strcmpi(gassType,"ang_farhang")
        phiFunc = @(X,e,phi,mu,alpha) alpha*phi + e*X;
        
    elseif strcmpi(gassType, "matthews_xie")
        phiFunc = @(X,e,phi,mu,alpha) e*X;
        
    else
        error(' invalid GASS type, must be of: \n {''benveniste'', ''ang_farhang'', ''matthews_xie''}')
    end
    
    % iterate over time
    for n=1:N
        % filter output n, y(n)
        y(n) = W(:, n)' * X(:, n);
        % prediction error n, e(n)
        e(n) = d(n) - y(n);
        % weights update rule
        W(:, n+1) = (1 - mu(n) * gamma) * W(:, n) + mu(n) * e(n) * X(:, n);
        % step-size update rule
        mu(n+1) = mu(n) + rho * e(n) * X(:, n)' * phi(:, n);
        % update phi
        phi(:, n+1) = phiFunc(X(:,n), e(n), phi(:,n), mu(n), alpha);
    end
    
    % discard first weight
    W = W(:, 2:end);
end