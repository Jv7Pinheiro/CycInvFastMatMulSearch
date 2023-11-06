function [G, fr, ft] = CI_Gradient_FunctionValue(CI_Mat, X, varargin)
%   CI_Gradient_FunctionValue returns the Gradient and the Function Value
%   in two seperate parts, fr is the real function value and fc is the 
%   corrected function through the targeted matrices.
%
%   G = CI_Gradient_FunctionValue(CI_Mat, X) returns the Gradient of X with
%   respect to the variables in CI_Mat implicitly.
%
%   G = CI_Gradient_FunctionValye(CI_Mat, X, varargin) specifies additional
%   parameters for the method. Specifically...
%
%   'beta' - correction parameter {0}
%   'TG_Mat' - Target Matrices {0 0 0 0}
%
%   Unless specified by input, the TG_Mat and lambda parameters are set to 0
%   so as to return the 'vanilla' Gradient and Function Value.
%
%   [G, fr] = CI_Gradient_FunctionValue(CI_Mat, X, varargin) also returns
%   the original Function Value.
%
%   [G, fr, ft] = CI_Gradient_FunctionValue(CI_Mat, X, varargin) also returns
%   the Function Value for targeted matrices part of the equation.

    %% Error Checking
    if (nargin < 2)
        error('Error: invalid input arguments');
    elseif numel(CI_Mat) ~= 4
        error('Error: Input cell array must contain exactly four elements.');
    elseif ndims(X) ~= 3
        error('Error: Tensor "X" is not 3-dimensional');
    end

    %% Set up parameters
    params = inputParser;
    params.addParameter('TG_Mat', {0 0 0 0},@(x) (iscell(x) || ismember(x,{'random','rand','randn','nvecs','zeros'})));
    params.addParameter('beta', [0 0 0 0], @(x) isvector(x) || isscalar(x));
    params.parse(varargin{:});
    
    %% Copy from params object
    TG_Mat = params.Results.TG_Mat;
    beta = params.Results.beta;
   
    %% Get from cells
    % Get matrices in cell
    S = CI_Mat{1};
    U = CI_Mat{2};
    V = CI_Mat{3};
    W = CI_Mat{4};

    % Get Sizes of tensor and chosen decompositions
    n = size(S, 1);
    Rs = size(S, 2);
    Rc = size(U, 2);

    if (Rs == 0 && Rc == 0)
        error("Rs and Rc cannot both be 0");
    end

    % Get target matrices in cell
    St = TG_Mat{1};
    Ut = TG_Mat{2};
    Vt = TG_Mat{3};
    Wt = TG_Mat{4};
    
    %% Precomputations
    % Transpose Precomputations
    if Rs ~= 0
        % S operations
        SS = S'*S;
    end

    if Rc ~= 0
        % U operations
        UU = U'*U;
        UV = U'*V;
        UW = U'*W;
        % V operations
        VU = V'*U;
        VV = V'*V;
        VW = V'*W;
        % W operations
        WU = W'*U;
        WV = W'*V;
        WW = W'*W;
    end
    
    if ~(Rs == 0 || Rc == 0)
        SU = S'*U;
        SV = S'*V;
        SW = S'*W;

        US = U'*S;
        VS = V'*S;
        WS = W'*S;
    end
    %% Get Gradients

    if Rs ~= 0 && Rc ~= 0
        dS = 3*(S*(SS.*SS) + U*(VS.*WS) + V*(US.*WS) + W*(US.*VS)) - (mttkrp(X, {S S S}, 1) + mttkrp(X, {S S S}, 2) + mttkrp(X, {S S S}, 3)) + beta(1)*(S - St);
        dU = 3*(S*(SV.*SW) + U*(VV.*WW) + V*(WV.*UW) + W*(UV.*VW)) - (mttkrp(X, {U W V}, 1) + mttkrp(X, {V U W}, 2) + mttkrp(X, {W V U}, 3)) + beta(2)*(U - Ut);
        dV = 3*(S*(SU.*SW) + U*(WU.*VW) + V*(UU.*WW) + W*(VU.*UW)) - (mttkrp(X, {V U W}, 1) + mttkrp(X, {W V U}, 2) + mttkrp(X, {U W V}, 3)) + beta(3)*(V - Vt);
        dW = 3*(S*(SU.*SV) + U*(VU.*WV) + V*(WU.*UV) + W*(UU.*VV)) - (mttkrp(X, {W V U}, 1) + mttkrp(X, {U W V}, 2) + mttkrp(X, {V U W}, 3)) + beta(4)*(W - Wt);

        G = [dS(:); dU(:); dV(:); dW(:)];
    else

        if (Rs == 0)
            % Compute everything that has U, V, W
            dU = 3*(U*(VV.*WW) + V*(WV.*UW) + W*(UV.*VW)) - (mttkrp(X, {U W V}, 1) + mttkrp(X, {V U W}, 2) + mttkrp(X, {W V U}, 3)) + beta(2)*(U - Ut);
            dV = 3*(U*(WU.*VW) + V*(UU.*WW) + W*(VU.*UW)) - (mttkrp(X, {V U W}, 1) + mttkrp(X, {W V U}, 2) + mttkrp(X, {U W V}, 3)) + beta(3)*(V - Vt);
            dW = 3*(U*(VU.*WV) + V*(WU.*UV) + W*(UU.*VV)) - (mttkrp(X, {W V U}, 1) + mttkrp(X, {U W V}, 2) + mttkrp(X, {V U W}, 3)) + beta(4)*(W - Wt);
            G = [dU(:); dV(:); dW(:)];
        elseif (Rc == 0)
            % Compute everything that has S
            dS = 3*(S*(SS.*SS)) - (mttkrp(X, {S S S}, 1) + mttkrp(X, {S S S}, 2) + mttkrp(X, {S S S}, 3)) + beta(1)*(S - St);

            G = dS(:);
        end
    end

    %% Get Function Value
    f1 = norm(X)^2;
    f2_1 = 0;
    f2_2 = 0;
    f2_3 = 0;
    f2_4 = 0;

    f3_1 = 0;
    f3_2 = 0;
    f3_3 = 0;
    f3_4 = 0;

    f4_1 = 0;
    f4_2 = 0;
    f4_3 = 0;
    f4_4 = 0;

    if Rs ~= 0
        f2_1 = dot(reshape(mttkrp(X, {S S S}, 1), [], 1), S(:));
        f3_1 = sum(reshape(SS.*SS.*SS, [], 1));
        f4_1 = norm(S-St)^2;
    end
    if Rc ~= 0
        f2_2 = dot(reshape(mttkrp(X, {U W V}, 1), [], 1), U(:));
        f2_3 = dot(reshape(mttkrp(X, {V U W}, 1), [], 1), V(:));
        f2_4 = dot(reshape(mttkrp(X, {W V U}, 1), [], 1), W(:));

        f3_3 = 3*sum(reshape(VV.*WW.*UU, [], 1));
        f3_4 = 6*sum(reshape(WV.*UW.*VU, [], 1));

        f4_2 = norm(U-Ut)^2;
        f4_3 = norm(V-Vt)^2;
        f4_4 = norm(W-Wt)^2;
    end
    if ~(Rs == 0 || Rc == 0)
        f3_2 = 6*sum(reshape(SV.*SW.*SU, [], 1));
    end

    f2 = f2_1 + f2_2 + f2_3 + f2_4;
    f3 = f3_1 + f3_2 + f3_3 + f3_4;
    
    % f's for target matrices function value
    ft = [f4_1 f4_2 f4_3 f4_4];
    
    % combine f_(1-3) 
    fr = 0.5*f1 - f2 + 0.5*f3;
end