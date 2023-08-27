%% Gradient
addpath("C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\tensor_toolbox-v3.5")
% Defining the model using symbolic vars
n = 4;
rank_s = 1;
rank_uvw = 2;
true_rank = rank_s + 3*rank_uvw;

Beta = sym('Beta', [1, 1], 'int');

S = sym('S', [n, rank_s], 'int');
U = sym('U',[n, rank_uvw], 'int');
V = sym('V',[n, rank_uvw], 'int');
W = sym('W',[n, rank_uvw], 'int');

St = sym('St', [n, rank_s], 'int');
Ut = sym('Ut',[n, rank_uvw], 'int');
Vt = sym('Vt',[n, rank_uvw], 'int');
Wt = sym('Wt',[n, rank_uvw], 'int');

X = sym('X',[n, n, n], 'int');

A = [S U V W];
B = [S W U V];
C = [S V W U];

% sum of outer products column vectors
X_hat = zeros(n, n, n);

for i=1:true_rank
    X_hat = X_hat + reshape(kron(kron(C(:,i), B(:,i)), A(:,i)), n, n, n);
end

phi = [X(:)-X_hat(:); Beta*(S(:)-St(:)); Beta*(U(:)-Ut(:)); Beta*(V(:)-Vt(:)); Beta*(W(:)-Wt(:))]; %phi
r = sum((phi .*phi)/2, "all");

%Compare Gradients
for i=1:1
    % randomize values to test gradient differences
    S_vals = randi(10, n, rank_s);
    U_vals = randi(10, n, rank_uvw);
    V_vals = randi(10, n, rank_uvw);
    W_vals = randi(10, n, rank_uvw);

    St_vals = randi(10, n, rank_s);
    Ut_vals = randi(10, n, rank_uvw);
    Vt_vals = randi(10, n, rank_uvw);
    Wt_vals = randi(10, n, rank_uvw);

    X_vals = randi(10, n, n, n);

    % put into factor matrices
    FM = {S_vals, U_vals, V_vals, W_vals};
    FtM = {St_vals, Ut_vals, Vt_vals, Wt_vals};

    %Just pass value matrices, and compute four gradients

    %compute gradient 1 from function
    [gradient1, f] = CI_Tilda_Gradient_FunctionValue(FM, FtM, tensor(X_vals), 9);

    % compute gradient 2 based on sym function and subbing variables
    gradient_sym = gradient(r, [S(:); U(:); V(:); W(:)]);
    gradient2 = subs(gradient_sym, S, S_vals);
    gradient2 = subs(gradient2, St, St_vals);
    gradient2 = subs(gradient2, U, U_vals);
    gradient2 = subs(gradient2, Ut, Ut_vals);
    gradient2 = subs(gradient2, V, V_vals);
    gradient2 = subs(gradient2, Vt, Vt_vals);
    gradient2 = subs(gradient2, W, W_vals);
    gradient2 = subs(gradient2, Wt, Wt_vals);
    gradient2 = subs(gradient2, X, X_vals);
    gradient2 = subs(gradient2, Beta, 3);
    gradient2 = cast(gradient2, "double");

    S_symval = reshape(gradient2(1:n*rank_s),n,rank_s);
    U_symval = reshape(gradient2(n*rank_s+1:n*rank_s+n*rank_uvw),n,rank_uvw);
    V_symval = reshape(gradient2(n*rank_s+1 + n*rank_uvw:n*rank_s+2*n*rank_uvw),n,rank_uvw);
    W_symval = reshape(gradient2(n*rank_s+1 +2*n*rank_uvw:n*rank_s+3*n*rank_uvw),n,rank_uvw);


    % message if gradients are ever differnt
    assert(isequal(gradient1, gradient2), "Gradients are different\n")
end
fprintf("Gradients are the same\n")


%% Jacobian
%Compare Jacobians
for i=1:100
    % randomize values to hessians
    S_vals = randi(10, n, rank_s);
    U_vals = randi(10, n, rank_uvw);
    V_vals = randi(10, n, rank_uvw);
    W_vals = randi(10, n, rank_uvw);
    X_vals = randi(10, n, n, n);

    % put into factor matrices
    FM = {S_vals, U_vals, V_vals, W_vals};

    %compute jacobian1 from our minres function
    jacobian1 = CI_Tilda_Jacobian(FM, 9);


    % compute jacobian 2 based on sym functions and subbing variables
    jacobian_sym = jacobian(phi(:), [S(:); U(:); V(:); W(:)]);

    jacobian2 = subs(jacobian_sym, S, S_vals);
    jacobian2 = subs(jacobian2, U, U_vals);
    jacobian2 = subs(jacobian2, V, V_vals);
    jacobian2 = subs(jacobian2, W, W_vals);
    jacobian2 = subs(jacobian2, Beta, -3);
    jacobian2 = subs(jacobian2, X, X_vals);
    jacobian2 = cast(jacobian2, "double");
    
    diff = jacobian1 + jacobian2;
    % message if jacobians are ever differnt
    assert(isequal(jacobian1, -jacobian2), "Jacobians are different")
end
fprintf("Jacobians are the same\n")