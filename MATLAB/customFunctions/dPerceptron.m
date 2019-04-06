function [y, e, W] = dPerceptron(X, d, mu, varargin)
% dPerceptron	Dynamic Perceptron, based on the LMS filter.
% Input: 
%       - X: Design matrix, [M N]
%       - d: Target vector, [1 N]
%       - mu: Step size, numeric
%       varargin:                                  variable input arguments
%       - varargin{1}   = gamma: leakage coefficient, scalar (default: 0)
%       - varargin{>1}  = options: bias, amplitude, pre-trained weights struct    (default: 0,1,[])
%       - varargin{any} = activator: activator, scalar       (default: @tanh)
% Output: 
%       * y: Filter output,     [1 N]
%       * e: Prediction error,  d-y
%       * W: Filter weights,    [M N]
% Usage: 
%   [y, e, w] = dPerceptron(X, d, mu, gamma) train LMS filter on Xd data.
%   [y, e, w] = dPerceptron(X, d, mu, @tanh) train dPerceptron on Xd 
%               data, with a tanh activator function and 0 bias, unscaled amplitude.
%   [y, e, w] = dPerceptron(X, d, mu, @tanh, opts) train dPerceptron on Xd 
%               using options specified in opts
%              *if pre-trained weights are specified that are the same
%               length as the design matrix/target vector - i.e. N, the
%               weights will not update

    gamma = 0; activator = @tanh; usePerceptron = false;
    opts.bias = 0; opts.ampl = 1; opts.W = [];  
    optsTemp=opts;
    updateWeights = true;

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
    
    % Check gamma, and usePerceptron
    if ~isempty(varargin)
        for ii=1:length(varargin)
            if isa(varargin{ii},'numeric') && ii==1 % gamma must be given first
                gamma = varargin{ii};
                if (gamma>1 && gamma<0)
                    warning('Gamma provided, %d, not between 0 and 1. \n Setting gamma to zero ', gamma)
                end
            elseif isstruct(varargin{ii})  % bias is given after gamma/activator function
                optsTemp = varargin{ii};
                usePerceptron = true;
            elseif isa(varargin{ii},'function_handle')      % activator function given anytime
                activator = varargin{ii};  % enable augmented CLMS - i.e. ACLMS
                usePerceptron = true;
            end
        end
    end
    
    optsFields = fieldnames(opts);
    for ii=find( isfield(optsTemp,optsFields) ==1 )'
        opts.(optsFields{ii}) = optsTemp.(optsFields{ii});
    end
    
    if opts.bias ~= 0
       X = [ones( 1, size(X,2) )*opts.bias ; X]; 
    end
    
    % sizes
    [M, N] = size(X);
    % Filter Output: pre-allocate for speed
    y = zeros(size(d));
    % Prediction Error: pre-allocate for speed
    e = zeros(size(d));
    % LMS Filter Weights: pre-allocate for speed
    W = zeros(M, N);        % use 0 weights to start with
    if ~isempty(opts.W)
        W(:,1:size(opts.W,2)) = opts.W;    % set initial weights to pre-trained
        if size(opts.W,2) == N
            updateWeights = false;         % if all weights set, don't update
            warning('weight update disbaled, fully-defined weights provided')
        end
    end
    
    % Iterate over the discrete time samples
    for n=1:N
        % Filter output n, y(n)
        if usePerceptron
            y(n) = opts.ampl *activator(  W(:,n)' *X(:,n) ) ;
        else
            y(n) = W(:,n)' *X(:,n);
        end
        % Prediction error n, e(n)
        e(n) = d(n) - y(n);
        % Weights update rule
        if updateWeights
            W(:,n+1) = (1 - mu*gamma) * W(:,n) + mu*e(n)*X(:, n);
        end
    end
    
    % Discard first weight
    W = W(:, 2:end);
    
    % Check if exploded
    if find(isnan(y)==1,1)
        warning('unstable mu provided, output reached NaN')
    end
end