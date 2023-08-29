function Q = CI_sparsify(P, threshold)
    
    %Q = P;
    dim = sqrt(size(P{1},1));
    Rs = size(P{1},2);
    Rc = size(P{2},2);
    
    %tmp
    T = matmul_tensor(dim,dim,dim);

    % check input error
    U = Convert_CImat_to_fac(P);
    inerr = norm(T-full(ktensor(U)))^2
    
    % count input nnz
    innnz = sum(cellfun(@(x) nnz(abs(x) > threshold), P))

    % try symmetric components
    Qs = cell(Rs,1);
    nnzs = zeros(Rs,1);
    for i = 1:Rs
        [Qs{i},~] = schur(reshape(P{1}(:,i),dim,dim));
        Q = cellfun(@(x) kron(Qs{i}',Qs{i}') * x, P, 'UniformOutput', false);
        nnzs(i) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q)); 
    end

    % try cyclic components
    Qc = cell(Rc,3);
    nnzc = zeros(Rc,3);
    for i = 1:Rc
        for j = 1:3
            [Qc{i,j},~] = schur(reshape(P{j+1}(:,i),dim,dim));
            Q = cellfun(@(x) kron(Qc{i,j}',Qc{i,j}') * x, P, 'UniformOutput', false);
            nnzc(i,j) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q)); 
        end
    end
    nnzs;
    nnzc;

    % choose sparsest solution
    [min_nnzs,i_nnzs] = min(nnzs);
    min_nnzc = min(nnzc(:));
    if min_nnzs < min_nnzc
        Q = cellfun(@(x) kron(Qs{i_nnzs}',Qs{i_nnzs}') * x, P, 'UniformOutput', false);
    else
        [i_nnzc,j_nnzc] = find(nnzc == min_nnzc,1);
        Q = cellfun(@(x) kron(Qc{i_nnzc,j_nnzc}', Qc{i_nnzc,j_nnzc}') * x, P, 'UniformOutput', false);
    end

    % check output error
    Unew = Convert_CImat_to_fac(Q);
    outerr = norm(T-(ktensor(Unew)))^2

    % count output nnz
    outnnz = sum(cellfun(@(x) nnz(abs(x) > threshold), Q))

end