function [K, P0, output] = ci_cp_dgn(Z, Rs, Rc, varargin)
%   CI_CP_OPT Fits a CP model to a cyclic invariant tensor via Damped Gauss-Newton optimization.
%
%   K = CI_CP_DGN(X, Rs, Rc) fits a cyclic invariant CANDECOMP/PARAFAC (CP) model
%   to the tensor X of chosen Rs and Rc length. The result K is a ktensor. The function being
%   optimized is F(K) = 1/2 || Z - K ||^2.
%
%   K = CI_CP_DGN(X, Rs, Rc,'param',value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'tol' - Tolerance on difference in relative error {1e-4}
%   'maxiters' - Maximum number of iterations {50}
%   'lambda' - Damping factor {1e-4}
%   'minstepsize' - Minimum stepsize for backtracking {1e-4}
%   'printitn' - Print fit every n iterations; 0 for no printing {1}
%   'targets' - Target matrices
%   'init' - Initialization for factor matrices (default: 'randn'). This
%   can be a cell array with the initial matrices, a ktensor, or one of the
%   following strings:
%      'randn'  Randomly generated via randn function
%      'rand'   Randomly generated via rand function
%      'zeros'  All zeros
%      'nvecs'  Selected as leading left singular vectors of X(n)
%
%   [K, U0] = CI_CP_DGN(...) also returns the initial guess.
%
%   [K, U0, OUT] = CI_CP_DGN(...) also returns a structure with information
%   on final relative error and related quantities.
    
    if (nargin < 3)
        error('Error: invalid input arguments');
    end

    %% Set parameters
    params = inputParser;
    params.addParameter('init', 'randn', @(x) (iscell(x) || isa(x, 'ktensor') || ismember(x,{'random','rand','randn','nvecs','zeros'})));
    params.addParameter('tol',1e-4,@isscalar);
    params.addParameter('maxiters',50,@(x) isscalar(x) & x > 0);
    params.addParameter('lambda',1e-4,@(x) isscalar(x) & x >= 0);
    params.addParameter('minstepsize',1e-4,@(x) isscalar(x) & x > 0);
    params.addParameter('printitn',1,@isscalar);
    params.addParameter('targets', {0 0 0 0}, @(x) iscell(x));
    params.addParameter('beta', [0 0 0 0], @(x) isvector(x) || isscalar(x));
    params.addParameter('cg_tol',1e-4,@isscalar);
    params.addParameter('cg_maxiters',20,@isscalar);
    params.parse(varargin{:});

    %% Copy from params object
    init = params.Results.init;
    tol = params.Results.tol;
    maxiters = params.Results.maxiters;
    lambda = params.Results.lambda;
    minstepsize = params.Results.minstepsize;
    printitn = params.Results.printitn;
    targets = params.Results.targets;
    beta = params.Results.beta;
    cg_tol = params.Results.cg_tol;
    cg_maxiters = params.Results.cg_maxiters;

    %% Initialization
    sz = size(Z, 1);

    % If there is a init parameter, then that becomres the initial guess,
    % else, a random first guess is made. Setting and saving the seed will
    % make answer reproduceble. 
    if iscell(init)
        if (size(init, 1) == 4)
            P0 = init;
        elseif(size(init, 1) == 1)
            P0 = init';
        end
    else
        P0 = cell(4,1);
        for n=1:4
            if (n==1)
                P0{n} = matrandnorm(feval(init, sz, Rs));
            else
                P0{n} = matrandnorm(feval(init, sz, Rc));
            end
        end
    end
    
    % If beta is a scalar, then a vector of length four where each entry is
    % the parameter beta. If vector is a length of two, then the first beta
    % will be applied to the S matrix and the second will be applied to U,
    % V, and W. If argument is a vector of length four, then nothing needs
    % to be done.
    if isscalar(beta)
        beta = repmat(beta, 4, 1)';
    elseif isvector(beta) 
        if length(beta) == 2
            beta = [beta(1) beta(2) beta(2) beta(2)];
        elseif length(beta) == 3 || length(beta) > 4
            error(['Beta must either be \n' ...
                'a single number that applies to all cyclic invariant matrices or\n' ...
                'a vector of length two where the first entry applies to S and the second entry applies to U, V, and W or \n' ...
                'a vector of length four where each entry corresponds to a cyclic invarient matrix'])
        end
    end


%% Fit CP using CPDGN

% The following output table contains:
% 'Iter' - Current cp_dgn iteration.
% 'f' - The combined function value of fr and ft. When beta = 0, f=fr.
% 'fr' - The original function.
% 'ft' - The targeted function value which includes beta and the target
% matrices.
% 'delta' - The change in the previous function value to new function value
% if this decreases then ci_cp_dgn breaks as progress is no longer being
% made
% 'stepsz' - 
% 'cgflag' - 
% 'cgiters' - 
% 'cgrelres' - 

if printitn > 0
    fprintf('\nCI_CP_DGN:\n');
    fprintf('Iter  |       f |      fr |      ft |  delta  |  stepsz |  cgflag | cgiters | cgrelres |\n');
    fprintf('----------------------------------------------------------------------------------------\n');
end

K = P0; % Set K to be initial guess
[G, fr, ft] = CI_Gradient_FunctionValue(P0, Z, 'beta', beta, 'TG_Mat', targets); % Compute function values and gradient of init

f = fr + beta(1)*ft(1) + beta(2)*ft(2) + beta(3)*ft(3) + beta(4)*ft(4); % Compute combined function value

for iter = 1:maxiters    
    % Compute search direction using CG
    [dvec,cg_flag,cg_relres,cg_iter] = pcg(@(x) Apply_Approx_Hessian(K, x, 'beta', beta, 'lambda', lambda), -G, cg_tol, cg_maxiters);

    D = CP_Vec_to_CI(dvec, K);
    if (Rs == 0)
        D{1} = 0;
    end
    % perform backtracking line search to ensure function value decrease
    alpha = 1;
    Kprev = K;
    fold = f;
    
    while alpha >= minstepsize

        % take Gauss-Newton step
        K = cellfun(@(x,y) x + alpha * y, Kprev, D, 'UniformOutput', false);
                
        [G, fr, ft] = CI_Gradient_FunctionValue(K, Z, 'beta', beta, 'TG_Mat', targets); % Compute func value and gradient of init
        f = fr + beta(1)*ft(1) + beta(2)*ft(2) + beta(3)*ft(3) + beta(4)*ft(4); % Compute combined function value

        delta = fold - f;
            
        % break if function value has decreased
        if delta > 0
            break
        end               
            
        % for next iteration decrease step size
        alpha = alpha / 2;
    end

    % Check for convergence
    if (iter > 1) && (delta < tol)
        flag = 0;
    else
        flag = 1;
    end
    
    % Print information after each iteration
    if ((mod(iter, printitn) == 0) || ((printitn > 0) && (flag == 0)))
        fprintf('%3d   | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e |   %3d   | %7.1e  | \n', iter, f, fr, sum(ft), delta, alpha, cg_flag, cg_iter, cg_relres);
    end
    
    % Check for convergence
    if (flag == 0)
        break;
    end 
end


output.FcnVal = f; % Return combined function value
output.fr = fr; % Return normal function value
output.ft = sum(ft); % Return the targeted function value
output.ExitFlag  = flag; % Return wheter it converged or reached max iter
output.NumIter = iter; % Return number of iterations
