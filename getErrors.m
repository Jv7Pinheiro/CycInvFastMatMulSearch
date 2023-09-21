function [rnd, rel, abs] = getErrors(K)
    dim = sqrt(size(K{1},1));
    
    T = matmul_tensor(dim,dim,dim);
    U = Convert_CI_to_Mat(K);
    R = Convert_CI_to_Mat(cellfun(@(x) round(x), K, 'UniformOutput', false));
    
    abs = norm(T - full(ktensor(U)));
    rel = norm(T - full(ktensor(U)))/norm(T);
    rnd = norm(T - full(ktensor(R)));
end