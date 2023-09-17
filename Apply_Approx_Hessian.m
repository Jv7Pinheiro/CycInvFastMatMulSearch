function JJx = Apply_Approx_Hessian(CI_Mat, x, varargin)
%   Apply_Approx_Hessian applies the approximate hessian J'J implicitly.
%
%   JJx = Apply_Approx_Hessian(CI_Mat, x) implicitly applies the 
%   approximate hessian using the cyclic invariant matrices in CI_Mat to
%   vector x.
%
%   JJx = Apply_Approx_Hessian(CI_Mat, x, param=value) specifies additional
%   parameters for the method. Specifically...
%
%   'beta' - correction parameter {0}
%   'lambda' - Damping factor {0}
%   
%   Unless specified by input, the beta and lambda parameters are set to 0
%   so as to apply the 'vanilla' approximate Hessian.

    %% Error Checking
    if (nargin < 2)
        error('Error: invalid input arguments');
    end

    %% Set up parameters
    params = inputParser;
    params.addParameter('beta', [0 0 0 0], @(x) isvector(x) || isscalar(x));
    params.addParameter('lambda', 0, @(x) isscalar(x) & x >= 0);
    params.parse(varargin{:});

    %% Copy from params object
    beta = params.Results.beta;
    lambda = params.Results.lambda;
   
    %% Get Values
    % Get matrices in cell
    S = CI_Mat{1};
    U = CI_Mat{2};
    V = CI_Mat{3};
    W = CI_Mat{4};
    
    % Get Sizes of tensor and chosen decompositions
    n = size(S, 1);
    Rs = size(S, 2);
    Rc = size(U, 2);
    
    % Get subvectors
    A = tt_cp_vec_to_ci_mat(x, CI_Mat);
    Ms = A{1};
    Mu = A{2};
    Mv = A{3};
    Mw = A{4};
    
    %% Transpose Precomputations
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
    % Ms operations
    MsS = Ms'*S;
    MsU = Ms'*U;
    MsV = Ms'*V;
    MsW = Ms'*W;
    % Mu operations
    MuS = Mu'*S;
    MuU = Mu'*U;
    MuV = Mu'*V;
    MuW = Mu'*W;
    % Mv operations
    MvS = Mv'*S;
    MvU = Mv'*U;
    MvV = Mv'*V;
    MvW = Mv'*W;
    % Mw operations
    MwS = Mw'*S;
    MwU = Mw'*U;
    MwV = Mw'*V;
    MwW = Mw'*W;
        
    %% Apply Approx Hessian to Vector 
    % If beta is 0 (which is the default option for ci_cp_dgn) it will not
    % affect the operations in any way and it would just return the normal
    % applied approximate hessian
    
    % Js Operations
    OneA = 3*(Ms*(SS.*SS) + 2*(S*(MsS.*SS))) + beta(1)*Ms; %Js'*Js*vec(Ms) + Bs*Ms
    OneB = 3*(Mu*(VS.*WS) + V*(MuS.*WS) + W*(MuS.*VS)); %Js'*Ju*vec(Mu)
    OneC = 3*(Mv*(US.*WS) + U*(MvS.*WS) + W*(MvS.*US)); %Js'*Jv*vec(Mv)
    OneD = 3*(Mw*(US.*VS) + U*(MwS.*VS) + V*(MwS.*US)); %Jw'*Jw*vec(Mw)
    % Ju Operations
    TwoA = 3*(Ms*(SV.*SW) + S*((MsW.*SV) + (MsV.*SW))); %Ju'*Js*vec(Ms)
    TwoB = 3*(Mu*(VV.*WW) + V*(MuW.*WV) + W*(MuV.*VW)) + beta(2)*Mu; %Ju'*Ju*vec(Mu) + Bu*Mu
    TwoC = 3*(Mv*(WV.*UW) + W*(MvW.*UV) + U*(MvV.*WW)); %Ju'*Jv*vec(Mv)
    TwoD = 3*(Mw*(UV.*VW) + U*(MwW.*VV) + V*(MwV.*UW)); %Ju'*Jw*vec(Mw)
    % Jv Operations
    ThreeA = 3*(Ms*(SU.*SW) + S*((MsU.*SW) + (MsW.*SU))); %Jv'*Js*vec(Ms)
    ThreeB = 3*(Mu*(VW.*WU) + V*(MuU.*WW) + W*(MuW.*VU)); %Jv'*Ju*vec(Mu)
    ThreeC = 3*(Mv*(WW.*UU) + W*(MvU.*UW) + U*(MvW.*WU)) + beta(3)*Mv; %Jv'*Jv*vec(Mv) + Bv*Mv
    ThreeD = 3*(Mw*(UW.*VU) + U*(MwU.*VW) + V*(MwW.*UU)); %Jv'*Jw*vec(Mw)
    %Jw Operations
    FourA = 3*(Ms*(SU.*SV) + S*((MsU.*SV) + (MsV.*SU))); %Jw'*Js*vec(Ms)
    FourB = 3*(Mu*(VU.*WV) + V*(MuV.*WU) + W*(MuU.*VV)); %Jw'*Ju*vec(Mu) 
    FourC = 3*(Mv*(WU.*UV) + W*(MvV.*UU) + U*(MvU.*WV)); %Jw'*Jv*vec(Mv)
    FourD = 3*(Mw*(UU.*VV) + U*(MwV.*VU) + V*(MwU.*UV)) + beta(4)*Mw; %Jw'*Jw*vec(Mw) + Bw*Mw
    
    % Join them together
    One = OneA + OneB + OneC + OneD + lambda*Ms;
    Two = TwoA + TwoB + TwoC + TwoD + lambda*Mu;
    Three = ThreeA + ThreeB + ThreeC + ThreeD + lambda*Mv;
    Four = FourA + FourB + FourC + FourD + lambda*Mw;
    
    % Return Vectorized result
    JJx = [One(:); Two(:); Three(:); Four(:)];
end