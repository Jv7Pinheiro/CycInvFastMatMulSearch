function A = ci_tt_cp_vec_to_fac(x,Mat)
% CI_TT_CP_VEC_TO_FAC Converts a vector to a cell array of the cyclic invariant matrices
% that are used to create the factor matrixes.
%
%   A = TT_CP_VEC_TO_FAC(X,Mat) converts the vector X into a cell array
%   of cyclic invariant matrices consistent with the size of the original matrices in Mat.
%
%


%% Set-up
N = size(Mat{1}, 1);
P = numel(x);
Rs = size(Mat{1}, 2);
Rc = size(Mat{2}, 2);

%if ~isequal(P, N*Rs + 3*N*Rc)
%    error('Sizes dont match')
%end

%% Create A
A = cell(4,1);
for n = 1:4
    if (n == 1)
        A{n} = reshape(x(1:N*Rs), [], Rs);
    else
        idx1 = N*Rs + (n-2)*N*Rc + 1;
        idx2 = N*Rs + (n-1)*N*Rc;
        A{n} = reshape(x(idx1:idx2), [], Rc);
    end
end