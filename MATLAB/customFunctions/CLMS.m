function [y, e, H, G] = CLMS(X, d, mu, varargin)
% CLMS	Complex Least Mean Square (CLMS) adaptive filter.
%       - X: design matrix, size(X)=[M N]
%       - d: target vector, size(d)=[1 N]
%       - mu: step size, scalar
%       varargin:                                  variable input arguments
%       - gamma: leakage coefficient, scalar
%       - setAug: 'aug' - enables ACLMS instead of CLMS
%       * y: filter output, size(y)=[1 N]
%       * e: prediction error, d(n) - y(n)
%       * H: filter weights, size(W)=[M N]
%       * G: filter weights, size(W)=[M N]
%   [y, e, H] = CLMS(X, d, mu) train CLMS filter on Xd data.
%   [y, e, H, G] = CLMS(X, d, mu,'aug') train ACLMS filter on Xd data.

    % Default values
    setAug = false; gamma = 0;

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
    
    % check if leakage coefficient is scalar and if augmented flag is true
    if ~isempty(varargin)
        for ii=1:length(varargin)
            if isscalar(varargin{ii})
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
    % filter output: init
    y = complex(zeros(size(d)));
    % prediction error: init
    e = complex(zeros(size(d)));
    % CLMS filter weights: init
    H = complex(zeros(M, N));
    % ACLMS filter weights: init
    G = complex(zeros(M, N));
    
    % NOTE: Matlab {'} automatically does complex conjugate traspose / Hermitian
    %       Matalab{.'} would do a pure transpose
    
    % iterate over time
    for n=1:N
        % filter output n, y(n)
        if setAug   % explicit difference to minimise computation
            y(n) = H(:,n)' *X(:,n) + G(:,n)' *conj( X(:,n) );
        else
            y(n) = H(:,n)' *X(:,n);
        end
        % prediction error n, e(n)
        e(n) = d(n) - y(n);
        % weights update rules
        if setAug   % explicit difference to minimise computation
            G(:,n+1) = (1- gamma*mu) *G(:,n) + mu*conj( e(n) )*conj( X(:,n) );
            H(:,n+1) = (1- gamma*mu) *H(:,n) + mu*conj( e(n) )*X(:,n);
        else
            H(:,n+1) = (1- gamma*mu) *H(:,n) + mu*conj( e(n) )*X(:,n);
        end
    end
    
    % discard first weight
    H = H(:, 2:end);
    G = G(:, 2:end);
end