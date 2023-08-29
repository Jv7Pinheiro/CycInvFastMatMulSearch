function [K, P0, output] = ci_cp_dgn(Z, Rs, Rc, varargin)
%   CI_CP_OPT Fits a CP model to a cyclic invariant tensor via Damped Gauss-Newton optimization.
%
%   K = CP_DGN(X,R) fits an R-component CANDECOMP/PARAFAC (CP) model
%   to the tensor X. The result K is a ktensor. The function being
%   optimized is F(K) = 1/2 || Z - K ||^2.
%
%   K = CP_DGN(X,R,'param',value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'tol' - Tolerance on difference in relative error {1e-4}
%   'maxiters' - Maximum number of iterations {50}
%   'lambda' - Damping factor {1e-4}
%   'minstepsize' - Minimum stepsize for backtracking {1e-4}
%   'printitn' - Print fit every n iterations; 0 for no printing {1}
%   'thresh' - Threshold for how low we want numbers to stay same.
%   'beta' - Something else
%   'init' - Initialization for factor matrices (default: 'randn'). This
%   can be a cell array with the initial matrices, a ktensor, or one of the
%   following strings:
%      'randn'  Randomly generated via randn function
%      'rand'   Randomly generated via rand function
%      'zeros'  All zeros
%      'nvecs'  Selected as leading left singular vectors of X(n)
%
%   [K, U0] = CP_DGN(...) also returns the initial guess.
%
%   [K, U0, OUT] = CP_DGN(...) also returns a structure with information
%   on final relative error and related quantities.
%
%   REFERENCE: (To be added)
%
%
%   See also TENSOR, SPTENSOR, KTENSOR.
%
%   MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation. (To be updated)


    %% Error checking
    % if ~isa(Z,'tensor') && ~isa(Z,'sptensor')
    %     error('Z must be a tensor or a sptensor');
    % end
    % 
    if (nargin < 3)
        error('Error: invalid input arguments');
    end

    %% Set parameters
    params = inputParser;
    params.addParameter('init', 'randn', @(x) (iscell(x) || isa(x, 'ktensor') || ismember(x,{'random','rand','randn','nvecs','zeros'})));
    params.addParameter('tol',1e-4,@isscalar);
    params.addParameter('maxiters',50,@(x) isscalar(x) & x > 0);
    params.addParameter('lambda',1e-4,@(x) isscalar(x) & x > 0);
    params.addParameter('minstepsize',1e-4,@(x) isscalar(x) & x > 0);
    params.addParameter('printitn',1,@isscalar);
    params.addParameter('thresh',0.3, @(x) isscalar(x) & x > 0);
    params.addParameter('beta', 0.1, @(x) isscalar(x) & x > 0);
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
    thresh = params.Results.thresh;
    beta = params.Results.beta;
    cg_tol = params.Results.cg_tol;
    cg_maxiters = params.Results.cg_maxiters;

    %% Initialization
    sz = size(Z, 1);
    if iscell(init)
        if (size(init, 1) == 4)
            P0 = init;
        elseif(size(init, 1) == 1)
            P0 = init';
        end
        TildaFunctions = 1;
    else
        TildaFunctions = 0;
        P0 = cell(4,1);
        for n=1:4
            if (n==1)
                P0{n} = matrandnorm(feval(init, sz, Rs));
            else
                P0{n} = matrandnorm(feval(init, sz, Rc));
            end
        end
    end

%% Fit CP using CPDGN
if printitn > 0
    fprintf('\nCI_CP_DGN:\n');
    fprintf('Iter  |       f |  relerr |  abserr |  delta  |  stepsz |  cgflag | cgiters | cgrelres |\n');
    fprintf('----------------------------------------------------------------------------------------\n');
end

K = P0;
normZsqr = norm(Z)^2; % precompute norm
if (TildaFunctions == 1)
    Kt = cellfun(@(x) SetTargets(x, thresh), K, 'UniformOutput', false)';
    [G,f] = CI_Tilda_Gradient_FunctionValue(P0, Kt, Z, beta); % compute func value and gradient of init
else
    [G,f] = CI_Gradient_FunctionValue(P0, Z); % compute func value and gradient of random
end

relerr = sqrt(2*f/normZsqr);

for iter = 1:maxiters    
    % compute search direction using CG
    % if there are initial values, we use the tilda functions
    if (TildaFunctions == 1)
        [dvec,cg_flag,cg_relres,cg_iter] = pcg(@(x) Apply_Approx_Tilda_Hessian(K, x, beta, lambda), -G, cg_tol, cg_maxiters);
    else
        [dvec,cg_flag,cg_relres,cg_iter] = pcg(@(x) Apply_Approx_Hessian(K, x, lambda), -G, cg_tol, cg_maxiters);
    end
    
    D = ci_tt_cp_vec_to_fac(dvec, K);
    
    % perform backtracking line search to ensure function value decrease
    alpha = 1;
    Kprev = K;
    relerrold = relerr;

    while alpha >= minstepsize

        % take Gauss-Newton step
        K = cellfun(@(x,y) x + alpha * y, Kprev, D, 'UniformOutput', false);
                
        % compute function value and gradient
        if (TildaFunctions == 1)
            [G,f] = CI_Tilda_Gradient_FunctionValue(K, Kt, Z, beta); % compute func value and gradient of init
        else
            [G,f] = CI_Gradient_FunctionValue(K, Z); % compute func value and gradient of random
        end

        relerr = sqrt(2*f/normZsqr);
        delta = relerrold - relerr;
            
        % break if function value has decreased
        if delta > 0
            break
        end               
            
        % for next iteration decrease step size
        alpha = alpha / 2;
    end

    
    abserr = norm(Z - ktensor(Convert_CImat_to_fac(K)))^2;
    % Check for convergence
    if (iter > 1) && (delta < tol)
        flag = 0;
    else
        flag = 1;
    end
    
    if (mod(iter, printitn) == 0) || ((printitn > 0) && (flag == 0))
        fprintf('%3d   | %7.1e | %7.1e | %.5f | %7.1e | %7.1e | %7.1e |   %3d   | %7.1e  | \n', iter, f, relerr, abserr, delta, alpha, cg_flag, cg_iter, cg_relres);
    end
    
    % Check for convergence
    if (flag == 0)
        break;
    end 
end

output.ExitFlag  = flag;
output.FcnVal = f;
output.RelErr = relerr;
output.AbsErr = abserr;
output.Fit = 100 * (1 - relerr);

end