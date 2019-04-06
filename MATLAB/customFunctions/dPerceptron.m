function [y, e, W] = dPerceptron(X, d, mu, varargin)
% dPerceptron	Dynamic Perceptron, based on the LMS filter.
%       - X: design matrix, size(X)=[M N]
%       - d: target vector, size(d)=[1 N]
%       - mu: step size, scalar
%       varargin:                                  variable input arguments
%       - varargin{1}   = gamma: leakage coefficient, scalar (default: 0)
%       - varargin{>1}  = options: bias, amplitude, pre-trained weights struct    (default: 0,1,[])
%       - varargin{any} = activator: activator, scalar       (default: @tanh)
%       * y: filter output, size(y)=[1 N]
%       * e: prediction error, d(n) - y(n)
%       * W: filter weights, size(W)=[M N]
%   [y, e, w] = dPerceptron(X, d, mu, gamma) train LMS filter on Xd data.
%   [y, e, w] = dPerceptron(X, d, mu, @tanh) train dPerceptron on Xd 
%               data, with a tanh activator function and 0 bias, unscaled amplitude.
%              *if pre-trained weights are specified that are the same
%               length as the design matrix/target vector - i.e. N, the
%               weights will not update

    gamma = 0; activator = @tanh; usePerceptron = false;
    opts.bias = 0; opts.ampl = 1; opts.W = [];  
    optsTemp=opts;
    updateWeights = true;

    % check if design matrix is 2D
    if ~ismatrix(X)
        error("design matrix must be 2D");
    end
    
    % check if target vector is 1D
    if ~isvector(d)
        error("target vector must be 1D");
    end
    
    % check X-d size consistency
    if size(X, 2) ~= size(d, 2)
        if size(X, 2) == size(d.', 2)
            d = d.';    % Using MATLAB {.'} operator to prevent conjugate transpose of complex data
            warning('auto-transposing design matrix data')
        else
            error("design matrix and target vector sizes are inconsistent");
        end
    end
    
    % check if step-size is scalar
    if ~isa(mu,'numeric')
        error("step-size parameter must be scalar");
    end
    
    % check if leakage coefficient is scalar
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
    % filter output: init
    y = zeros(size(d));
    % prediction error: init
    e = zeros(size(d));
    % LMS filter weights: init
    W = zeros(M, N);        % use 0 weights to start with
    if ~isempty(opts.W)
        W(:,1:size(opts.W,2)) = opts.W;    % set initial weights to pre-trained
        if size(opts.W,2) == N
            updateWeights = false;         % if all weights set, don't update
            warning('weight update disbaled, fully-defined weights provided')
        end
    end
    
    % iterate over time
    for n=1:N
        % filter output n, y(n)
        if usePerceptron
            y(n) = opts.ampl *activator(  W(:,n)' *X(:,n) ) ;
        else
            y(n) = W(:,n)' *X(:,n);
        end
        % prediction error n, e(n)
        e(n) = d(n) - y(n);
        % weights update rule
        if updateWeights
            W(:,n+1) = (1 - mu*gamma) * W(:,n) + mu*e(n)*X(:, n);
        end
    end
    
    % discard first weight
    W = W(:, 2:end);
    
    if find(isnan(y)==1,1)
        warning('unstable mu provided, output reached NaN')
    end
end