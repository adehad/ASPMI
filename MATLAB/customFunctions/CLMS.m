function [y, e, H, G] = CLMS(X, d, mu, varargin)
% CLMS	Complex Least Mean Square (CLMS) adaptive filter.
% Input: 
%       - X: Design matrix, [M N]
%       - d: Target vector, [1 N]
%       - mu: Step size, numeric
%       varargin:                                  variable input arguments
%       - gamma: Leakage Coefficient, numeric
%       - setAug: 'aug' - enables ACLMS instead of CLMS
% Output: 
%       * y: Filter output,     [1 N]
%       * e: Prediction error,  d-y
%       * H: Filter weights,    [M N]
%       * G: Filter weights,    [M N]
% Usage: 
%   [y, e, H] = CLMS(X, d, mu) train CLMS filter on Xd data.
%   [y, e, H, G] = CLMS(X, d, mu,'aug') train ACLMS filter on Xd data.

    % Default values
    setAug = false; gamma = 0;

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
            error("Design matrix and target vector sizes are  incompatible, [M N] and [1 N] required");
        end
    end
    % Step-size is a numeric scalar
    if ~isa(mu,'numeric')
        error("Step-size parameter (mu) must be numeric");
    end
    
    % Check if leakage coefficient is numeric and if augmented flag is true
    if ~isempty(varargin)
        for ii=1:length(varargin)
            if isa(varargin{ii},'numeric')
                gamma = varargin{ii};
                if (gamma>1 && gamma<0)
                    warning('Gamma provided, %d, not between 0 and 1. \n Setting gamma to zero ', gamma)
                end
            elseif isstring(varargin{ii}) || ischar(varargin{ii})
                if strcmpi(varargin{ii},'aug')
                    setAug = true;  % enable augmented CLMS - i.e. ACLMS
                else
                    warning('Only string allowed is ''aug'' to activated ACLMS')
                end
            end
        end
    end
    
    % sizes
    [M, N] = size(X);
    % Filter Output: pre-allocate for speed
    y = complex(zeros(size(d)));
    % Prediction Error: pre-allocate for speed
    e = complex(zeros(size(d)));
    % CLMS Filter Weights: pre-allocate for speed
    H = complex(zeros(M, N));
    % ACLMS filter weights: pre-allocate for speed
    G = complex(zeros(M, N));
    
    % NOTE: Matlab {'} automatically does complex conjugate traspose / Hermitian
    %       Matalab{.'} would do a pure transpose
    
    % Iterate over the discrete time samples
    for n=1:N
        % Filter output n, y(n)
        if setAug   % explicit difference to minimise computation
            y(n) = H(:,n)' *X(:,n) + G(:,n)' *conj( X(:,n) );
        else
            y(n) = H(:,n)' *X(:,n);
        end
        % Prediction error n, e(n)
        e(n) = d(n) - y(n);
        % Weights update rules
        if setAug   % explicit difference to minimise computation
            G(:,n+1) = (1- gamma*mu) *G(:,n) + mu*conj( e(n) )*conj( X(:,n) );
            H(:,n+1) = (1- gamma*mu) *H(:,n) + mu*conj( e(n) )*X(:,n);
        else
            H(:,n+1) = (1- gamma*mu) *H(:,n) + mu*conj( e(n) )*X(:,n);
        end
    end
    
    % Discard first weight
    H = H(:, 2:end);
    G = G(:, 2:end);
    
    % Check Instability
    if find(isnan(y)==1,1)
        warning('unstable mu provided, output reached NaN')
    end
end