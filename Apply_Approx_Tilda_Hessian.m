function JJW = Apply_Approx_Tilda_Hessian(Mat, M, Beta, lambda)
    % Get matrices in cell
    S = Mat{1};
    U = Mat{2};
    V = Mat{3};
    W = Mat{4};
    
    n = size(S, 1);
    Rs = size(S, 2);
    Rc = size(U, 2);

    % Get subvectors
    A = ci_tt_cp_vec_to_fac(M,Mat);
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
    OneA = 3*(Ms*(SS.*SS) + 2*(S*(MsS.*SS))) + Beta*Ms; %Js'*Js*vec(Ms)
    OneB = 3*(Mu*(VS.*WS) + V*(MuS.*WS) + W*(MuS.*VS)); %Js'*Ju*vec(Mu)
    OneC = 3*(Mv*(US.*WS) + U*(MvS.*WS) + W*(MvS.*US)); %Js'*Jv*vec(Mv)
    OneD = 3*(Mw*(US.*VS) + U*(MwS.*VS) + V*(MwS.*US)); %Jw'*Jw*vec(Mw)
    % Ju Operations
    TwoA = 3*(Ms*(SV.*SW) + S*((MsW.*SV) + (MsV.*SW))); %Ju'*Js*vec(Ms)
    TwoB = 3*(Mu*(VV.*WW) + V*(MuW.*WV) + W*(MuV.*VW)) + Beta*Mu; %Ju'*Ju*vec(Mu)
    TwoC = 3*(Mv*(WV.*UW) + W*(MvW.*UV) + U*(MvV.*WW)); %Ju'*Jv*vec(Mv)
    TwoD = 3*(Mw*(UV.*VW) + U*(MwW.*VV) + V*(MwV.*UW)); %Ju'*Jw*vec(Mw)
    % Jv Operations
    ThreeA = 3*(Ms*(SU.*SW) + S*((MsU.*SW) + (MsW.*SU))); %Jv'*Js*vec(Ms)
    ThreeB = 3*(Mu*(VW.*WU) + V*(MuU.*WW) + W*(MuW.*VU)); %Jv'*Ju*vec(Mu)
    ThreeC = 3*(Mv*(WW.*UU) + W*(MvU.*UW) + U*(MvW.*WU)) + Beta*Mv; %Jv'*Jv*vec(Mv)
    ThreeD = 3*(Mw*(UW.*VU) + U*(MwU.*VW) + V*(MwW.*UU)); %Jv'*Jw*vec(Mw)
    %Jw Operations
    FourA = 3*(Ms*(SU.*SV) + S*((MsU.*SV) + (MsV.*SU))); %Jw'*Js*vec(Ms)
    FourB = 3*(Mu*(VU.*WV) + V*(MuV.*WU) + W*(MuU.*VV)); %Jw'*Ju*vec(Mu) 
    FourC = 3*(Mv*(WU.*UV) + W*(MvV.*UU) + U*(MvU.*WV)); %Jw'*Jv*vec(Mv)
    FourD = 3*(Mw*(UU.*VV) + U*(MwV.*VU) + V*(MwU.*UV)) + Beta*Mw; %Jw'*Jw*vec(Mw)
    
    % Join them together
    One = OneA + OneB + OneC + OneD + lambda*Ms;
    Two = TwoA + TwoB + TwoC + TwoD + lambda*Mu;
    Three = ThreeA + ThreeB + ThreeC + ThreeD + lambda*Mv;
    Four = FourA + FourB + FourC + FourD + lambda*Mw;
    

    % Return Vectorized result
    JJW = [One(:); Two(:); Three(:); Four(:);];
end