function JJW = Apply_Approx_Hessian(Mat, M, lambda)
    % Get matrices in cell
    S = Mat{1};
    U = Mat{2};
    V = Mat{3};
    W = Mat{4};
    
    % Get subvectors
    A = CyclicInvariant_tt_cp_vec_to_fac(M,Mat);
    Ms = A{1};
    Mu = A{2};
    Mv = A{3};
    Mw = A{4};

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
    % M operations
    % Z operations
    MsS = Ms'*S;
    MsU = Ms'*U;
    MsV = Ms'*V;
    MsW = Ms'*W;
    % X operations
    MuS = Mu'*S;
    MuU = Mu'*U;
    MuV = Mu'*V;
    MuW = Mu'*W;
    % Y operations
    MvS = Mv'*S;
    MvU = Mv'*U;
    MvV = Mv'*V;
    MvW = Mv'*W;
    % K operations
    MwS = Mw'*S;
    MwU = Mw'*U;
    MwV = Mw'*V;
    MwW = Mw'*W;
        
    % Apply Approx Hessian to Vector 
    % Js Operations
    OneA = 3*(Ms*(SS.*SS) + 2*(S*(MsS.*SS))); % Correct
    OneB = 3*(Mu*(VS.*WS) + V*(MuS.*WS) + W*(MuS.*VS)); % Correct
    OneC = 3*(Mv*(US.*WS) + U*(MvS.*WS) + W*(MvS.*US)); % Correct
    OneD = 3*(Mw*(US.*VS) + U*(MwS.*VS) + V*(MwS.*US)); % Correct
    % Ju Operations
    TwoA = 3*(Ms*(SV.*SW) + S*((MsW.*SV) + (MsV.*SW))); % Correct
    TwoB = 3*(Mu*(VV.*WW) + V*(MuW.*WV) + W*(MuV.*VW));
    TwoC = 3*(Mv*(WV.*UW) + W*(MvW.*UV) + U*(MvV.*WW));
    TwoD = 3*(Mw*(UV.*VW) + U*(MwW.*VV) + V*(MwV.*UW));
    % Jv Operations
    ThreeA = 3*(Ms*(SU.*SW) + S*((MsU.*SW) + (MsW.*SU)));
    ThreeB = 3*(Mu*(VW.*WU) + V*(MuU.*WW) + W*(MuW.*VU));
    ThreeC = 3*(Mv*(WW.*UU) + W*(MvU.*UW) + U*(MvW.*WU));
    ThreeD = 3*(Mw*(UW.*VU) + U*(MwU.*VW) + V*(MwW.*UU));
    %Jw Operations
    FourA = 3*(Ms*(SU.*SV) + S*((MsU.*SV) + (MsV.*SU))); % Correct
    FourB = 3*(Mu*(VU.*WV) + V*(MuV.*WU) + W*(MuU.*VV));
    FourC = 3*(Mv*(WU.*UV) + W*(MvV.*UU) + U*(MvU.*WV));
    FourD = 3*(Mw*(UU.*VV) + U*(MwV.*VU) + V*(MwU.*UV));
    
    % Join them together
    One = OneA + OneB + OneC + OneD + lambda*Ms;
    Two = TwoA + TwoB + TwoC + TwoD + lambda*Mu;
    Three = ThreeA + ThreeB + ThreeC + ThreeD + lambda*Mv;
    Four = FourA + FourB + FourC + FourD + lambda*Mw;
    % Return Vectorized result
    JJW = [One(:); Two(:); Three(:); Four(:)];
end