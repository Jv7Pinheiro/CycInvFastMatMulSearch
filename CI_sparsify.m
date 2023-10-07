function [Q, innnz, outnnz] = CI_sparsify(P, threshold, varargin)
    % Introduces zeros while maintaining the same error. Default type
    % uses schur decomposition.
    
    %% Set up parameters
    params = inputParser;
    params.addParameter('type', 'schur', @(x) isstring(x))
    params.addParameter('print', 0, @isscalar);
    params.parse(varargin{:});

    %% Copy from params object
    type = params.Results.type;
    print = params.Results.print;

    type = lower(type);
    if ~(isequal(type, 'eig') || isequal(type, 'schur'))
        error('function type must be either eig or schur');
    end

    %Q = P;
    dim = sqrt(size(P{1},1));
    Rs = size(P{1},2);
    Rc = size(P{2},2);
    
    %tmp
    T = matmul_tensor(dim,dim,dim);
    
    % check input error
    U = Convert_CI_to_Mat(P);
    inerr = norm(T-full(ktensor(U)))^2;
    
    % count input nnz
    innnz = sum(cellfun(@(x) nnz(abs(x) > threshold), P));

    if (print == 1)
        fprintf('\nInitial error is %7.1e, with %d zeros\n', inerr, innnz);
    end

    % try symmetric components
    Qs = cell(Rs,1);
    nnzs = zeros(Rs,1);
    if isequal(type, 'schur')
        for i = 1:Rs
            [Qs{i},~] = schur(reshape(P{1}(:,i),dim,dim));
            Q = cellfun(@(x) kron(Qs{i}',Qs{i}') * x, P, 'UniformOutput', false);
            nnzs(i) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q)); 
        end
    elseif (isequal(type, 'eig'))
        for i = 1:Rs
            PP = (reshape(P{1}(:, i),dim,dim)' + reshape(P{1}(:, i),dim,dim))/2;
            [Qs{i},~] = eig(reshape(PP,dim,dim));
            Q = cellfun(@(x) kron(Qs{i}',Qs{i}') * x, P, 'UniformOutput', false);
            nnzs(i) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q)); 
        end
    end

    % try cyclic components
    Qc = cell(Rc,3);
    nnzc = zeros(Rc,3);
    if isequal(type,'schur')
        for i = 1:Rc
            for j = 1:3
                [Qc{i,j},~] = schur(reshape(P{j+1}(:,i),dim,dim));
                Q = cellfun(@(x) kron(Qc{i,j}',Qc{i,j}') * x, P, 'UniformOutput', false);
                nnzc(i,j) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q));
            end
        end
    elseif isequal(type, 'eig')
        for i = 1:Rc
            for j = 1:3
                PP = (reshape(P{j+1}(:, i),dim,dim)' + reshape(P{j+1}(:, i),dim,dim))/2;
                [Qc{i,j},~] = eig(reshape(PP,dim,dim));
                Q = cellfun(@(x) kron(Qc{i,j}',Qc{i,j}') * x, P, 'UniformOutput', false);
                nnzc(i,j) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q));
            end
        end
    end

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
    Unew = Convert_CI_to_Mat(Q);
    outerr = norm(T-(ktensor(Unew)))^2;

    % count output nnz
    outnnz = sum(cellfun(@(x) nnz(abs(x) > threshold), Q));

    if (print == 1)
        fprintf('Output error is %7.1e, with %d zeros\n\n', outerr, outnnz);
    end

end