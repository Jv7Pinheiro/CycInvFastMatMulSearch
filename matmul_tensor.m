function T = matmul_tensor(m,n,p) 
% build m x n x p matmul tensor
% assumes col-major order of inputs and row-major order of output

    % init tensor
    T=tensor(zeros(m*n,n*p,m*p));
    
    % iterate over nonzero entries in tensor
    for i = 1:m
        for j = 1:n
            for k = 1:p
                T(i+(j-1)*m,j+(k-1)*n,k+(i-1)*p) = 1;
            end
        end
    end
    
end