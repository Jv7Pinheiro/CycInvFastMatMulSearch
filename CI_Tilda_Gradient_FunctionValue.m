function [df, f] = CI_Tilda_Gradient_FunctionValue(Mat, Tilda_Mat, X, Beta)
    %Check for correct input arguments
    if numel(Mat) ~= 4
        error('Input cell array must contain exactly four elements.');
    elseif ndims(X) ~= 3
        error('Tensor "X" is not 3-dimensional');
    end
    
    % Get matrices in cell
    S = Mat{1};
    U = Mat{2};
    V = Mat{3};
    W = Mat{4};

    % Get tilda matrices in cell
    St = Tilda_Mat{1};
    Ut = Tilda_Mat{2};
    Vt = Tilda_Mat{3};
    Wt = Tilda_Mat{4};
    
    %% Precomputations
    % Transpose Precomputations
    % S operations
    SS = S'*S;
    SU = S'*U;
    SV = S'*V;
    SW = S'*W;
    % U operations
    US = U'*S;
    UU = U'*U;
    UV = U'*V;
    UW = U'*W;
    % V operations
    VS = V'*S;
    VU = V'*U;
    VV = V'*V;
    VW = V'*W;
    % W operations
    WS = W'*S;
    WU = W'*U;
    WV = W'*V;
    WW = W'*W;
    
    %% Get Gradients
    dS = 3*(S*(SS.*SS) + U*(VS.*WS) + V*(US.*WS) + W*(US.*VS)) - (mttkrp(X, {S S S}, 1) + mttkrp(X, {S S S}, 2) + mttkrp(X, {S S S}, 3)) + Beta*(S - St);
    dU = 3*(S*(SV.*SW) + U*(VV.*WW) + V*(WV.*UW) + W*(UV.*VW)) - (mttkrp(X, {U W V}, 1) + mttkrp(X, {V U W}, 2) + mttkrp(X, {W V U}, 3)) + Beta*(U - Ut);
    dV = 3*(S*(SU.*SW) + U*(WU.*VW) + V*(UU.*WW) + W*(VU.*UW)) - (mttkrp(X, {V U W}, 1) + mttkrp(X, {W V U}, 2) + mttkrp(X, {U W V}, 3)) + Beta*(V - Vt);
    dW = 3*(S*(SU.*SV) + U*(VU.*WV) + V*(WU.*UV) + W*(UU.*VV)) - (mttkrp(X, {W V U}, 1) + mttkrp(X, {U W V}, 2) + mttkrp(X, {V U W}, 3)) + Beta*(W - Wt);

    df = [dS(:); dU(:); dV(:); dW(:)];

    %% Get Function Value
    f1 = norm(X)^2;

    f2_1 = dot(reshape(mttkrp(X, {S S S}, 1), [], 1), S(:));
    f2_2 = dot(reshape(mttkrp(X, {U W V}, 1), [], 1), U(:));
    f2_3 = dot(reshape(mttkrp(X, {V U W}, 1), [], 1), V(:));
    f2_4 = dot(reshape(mttkrp(X, {W V U}, 1), [], 1), W(:));
    f2 = f2_1 + f2_2 + f2_3 + f2_4;

    f3_1 = sum(reshape(SS.*SS.*SS, [], 1));
    f3_2 = 6*sum(reshape(SV.*SW.*SU, [], 1));
    f3_3 = 3*sum(reshape(VV.*WW.*UU, [], 1));
    f3_4 = 6*sum(reshape(WV.*UW.*VU, [], 1));
    f3 = f3_1 + f3_2 + f3_3 + f3_4;
    
    f4_1 = norm(S-St)^2;
    f4_2 = norm(U-Ut)^2;
    f4_3 = norm(V-Vt)^2;
    f4_4 = norm(W-Wt)^2;
    AddErr = f4_1 + f4_2 + f4_3 +f4_4
    f4 = Beta*AddErr
    
    ResErr = 0.5*f1 - f2 + 0.5*f3
    f = ResErr + f4;
end