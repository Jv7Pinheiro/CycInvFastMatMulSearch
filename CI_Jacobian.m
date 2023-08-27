function J = CI_Jacobian(Mat)
    %Check for correct input argument
    if numel(Mat) ~= 4
        error('*Input cell array must contain exactly four elements.');
    end

    S = Mat{1};
    U = Mat{2};
    V = Mat{3};
    W = Mat{4};

    n = size(S, 1);
    I = eye(n);

    % Precomputations
    SoS = khatrirao(S, S);
    VoW = khatrirao(V, W);
    WoV = khatrirao(W, V);
    WoU = khatrirao(W, U);
    UoW = khatrirao(U, W);
    UoV = khatrirao(U, V);
    VoU = khatrirao(V, U);

    % Permutation Matrices
    P2 = kron(eye(n),perfect_shuffle(n,n));
    P3 = perfect_shuffle(n^2,n);
    
    Js = kron(SoS, I) + P2'*kron(SoS, I) + P3'*kron(SoS, I);
    Ju = kron(VoW, I) + P2'*kron(WoV, I) + P3'*kron(VoW, I);
    Jv = kron(WoU, I) + P2'*kron(UoW, I) + P3'*kron(WoU, I);
    Jw = kron(UoV, I) + P2'*kron(VoU, I) + P3'*kron(UoV, I);

    J = [Js Ju Jv Jw];
end
