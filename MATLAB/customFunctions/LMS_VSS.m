function [y, e, W] = LMS_VSS(X, d, mu_0, gamma, rho, gassType, varargin)
% GASS	Gradient Adaptive Step-Size (GASS) with Least Mean Square (LMS) adaptive filter.
% Input: 
%       - X: Design matrix, [M N]
%       - d: Target vector, [1 N]
%       - mu_0: Initial step size, numeric
%       - gamma: leakage coefficient, numeric
%       - rho: learning rate, numeric
%       varargin:                                  variable input arguments
%       - gassType: algorithm name, string from
%               {'benveniste', 'ang_farhang', 'matthews_xie'}
%       - alpha: Ang & Farhang learning parameter, scalar
% Output: 
%       * y: Filter output,     [1 N]
%       * e: Prediction error,  d-y
%       * W: Filter weights,    [M N]
% Usage: 
%   [y, e, W] = LMS_VSS(X, d, mu_0, gamma, rho, gassType, alpha) train GASS filter on Xd data.
%               alpha only required for Ang & Farhang

    alpha = 0; % if not using the Ang & Farhang method, we don't mind if it is zero

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
    
    % Initial Step-size is a numeric scalar
    if ~isa(mu_0,'numeric')
        error("Step-size parameter (mu) must be numeric");
    end
    
    % Check if leakage coefficient is scalar
    if ~isscalar(gamma)
        error("Leakage coefficient parameter must be scalar");
    end
    
    % Check if learning rate is scalar
    if ~isa(rho,'numeric')
        error("Learning rate parameter must be scalar");
    end
    
    % Check if algorithm type is string
    if ~isstring(gassType)
        error("Algorithm type parameter must be string. e.g. ''benveniste'' OR ""benveniste"" ");
    end
    
    % Check if Ang & Farhang learning parameter is set and is scalar
    if ~isempty(varargin)
        alpha = varargin{1};
        if ~isa(alpha,'numeric')
            error("Ang & Farhang learning parameter (alpha) must be scalar");
        end
    else
        if strcmpi(gassType,"ang_farhang")
            error("Ang & Farhang learning parameter (alpha) must be specified");
        end
    end

    % sizes
    [M, N] = size(X);
    % Filter Output: pre-allocate for speed
    y = zeros(size(d));
    % Prediction Error: pre-allocate for speed
    e = zeros(size(d));
    % LMS VSS filter weights: pre-allocate for speed
    W = zeros(M, N+1);
    % Step-Size: pre-allocate for speed
    mu = zeros(1, N+1);
    mu(1) = mu_0;
    % phi term: pre-allocate for speed
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
    
    % Iterate over the discrete time samples
    for n=1:N
        % Filter output n, y(n)
        y(n) = W(:,n)' * X(:,n);
        % Prediction error n, e(n)
        e(n) = d(n) - y(n);
        % Weights update rule
        W(:,n+1) = (1 - mu(n)*gamma) * W(:,n) + mu(n) * e(n) * X(:,n);
        % Step-Size update rule
        mu(n+1) = mu(n) + rho * e(n) * X(:,n)' * phi(:,n);
        % Update phi
        phi(:,n+1) = phiFunc( X(:,n), e(n), phi(:,n), mu(n), alpha );
    end
    
    % Discard first weight
    W = W(:,2:end);
    
    
    % Check Instability
    if find(isnan(y)==1,1)
        warning('unstable mu provided, output reached NaN')
    end
end