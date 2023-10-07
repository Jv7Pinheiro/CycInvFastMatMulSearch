%% Set Up
MMT2 = matmul_tensor(2, 2, 2); 
% Rank 7 | 2 Decompositions
% (Rs = 4, Rc = 1), (Rs = 1, Rc = 2).

MMT3 = matmul_tensor(3, 3, 3);
% Rank 23 | 7 Decompositions
% (Rs = 20, Rc = 1), (Rs = 17, Rc = 2), (Rs = 14, Rc = 3),
% (Rs = 11, Rc = 4), (Rs = 8, Rc = 5), (Rs = 5, Rc = 6), 
% (Rs = 2, Rc = 7)
% Rank 22 | 7 Decompositions
% (Rs = 19, Rc = 1), (Rs = 16, Rc = 2), (Rs = 13, Rc = 3),
% (Rs = 10, Rc = 4), (Rs = 7, Rc = 5), (Rs = 4, Rc = 6), 
% (Rs = 1, Rc = 7)
% Rank 21 | 6 Decompositions
% (Rs = 18, Rc = 1), (Rs = 15, Rc = 2), (Rs = 12, Rc = 3),
% (Rs = 9, Rc = 4), (Rs = 6, Rc = 5), (Rs = 3, Rc = 6)
% Rank 20 | 6 Decompositions
% (Rs = 17, Rc = 1), (Rs = 14, Rc = 2), (Rs = 11, Rc = 3),
% (Rs = 8, Rc = 4), (Rs = 5, Rc = 5), (Rs = 2, Rc = 6)
% Rank 19 | 6 Decompositions
% (Rs = 16, Rc = 1), (Rs = 13, Rc = 2), (Rs = 10, Rc = 3),
% (Rs = 7, Rc = 4), (Rs = 4, Rc = 5), (Rs = 1, Rc = 6)

MMT4 = matmul_tensor(4, 4, 4); 
% Rank 49 | 16 Decompositions
% (Rs = 46, Rc = 1), (Rs = 43, Rc = 2), (Rs = 40, Rc = 3),
% (Rs = 37, Rc = 4), (Rs = 34, Rc = 5), (Rs = 31, Rc = 6), 
% (Rs = 28, Rc = 7), (Rs = 25, Rc = 8), (Rs = 22, Rc = 9),
% (Rs = 19, Rc = 10), (Rs = 16, Rc = 11), (Rs = 13, Rc = 12), 
% (Rs = 10, Rc = 13), (Rs = 7, Rc = 14), (Rs = 4, Rc = 15), 
% (Rs = 1, Rc = 16)

% Set which tensor to test
T = MMT2; % Decomposing Tensor

if isequal(T, MMT2)
    NumItr = 50;
    Tensor = 'MMT2';
    
    Rs = 1;
    Rc = 2;
    
    clear MMT3 MMT4;
elseif isequal(T, MMT3)
    NumItr = 5000;
    Tensor = 'MMT3';

    Rs = 2;
    Rc = 7;

    clear MMT2 MMT4;
elseif isequal(T, MMT4)
    NumItr = 10000;
    Tensor = 'MMT4';
    
    Rs = 1;
    Rc = 16;

    clear MMT2 MMT3;
end

% Thresholds used in CI_sparsify, roundWithThresholds, and Setting Targets
thresh = [0.1 0.2 0.3 0.4 0.5]';
t_sz = size(thresh, 1);

% This is how many times we will use a out_cell as in_cell
MaxOuterItr = 5;

%% Start Testing
Data = zeros(NumItr, t_sz, MaxOuterItr, 6);
fprintf('\nSearching %s solutions with Rs=%d, Rc=%d\n', Tensor, Rs, Rc);

% str = cell(1, NumItr);
% % Initialize each random number stream with a unique seed
% for i = 1:NumItr
%     str{i} = RandStream.create('mlfg6331_64', 'Seed', i);
% end

tic;
parfor i = 1:NumItr
    for j = 1:t_sz
        for k = 1:MaxOuterItr
            if k == 1
                rng(i, 'combRecursive');
                %seed = rng().Seed;
                % RandStream.setGlobalStream(str{i});
                % seed = RandStream.getGlobalStream.Seed
                % Run initial test with randomized values from seed
                K = ci_cp_dgn(T, Rs, Rc, 'printitn', 0, 'maxiters', 150, 'lambda', 1e-6, 'tol', 1e-8);
            else
                K = ci_cp_dgn(T, Rs, Rc, 'init', RSP_K, 'printitn', 0, 'maxiters', 150, 'lambda', 1e-6, 'tol', 1e-8);
            end
            [rnd_cp, ~, abs_cp] = getErrors(K);
            [SP_K, innz, outnz] = CI_sparsify(K, thresh(j));
            % Rounding with Threshold
            RSP_K = cellfun(@(x) roundWithThreshold(x, thresh(j)), SP_K, 'UniformOutput', false); 
            
            [rnd_rsp, ~, abs_rspp] = getErrors(RSP_K);
            Data(i, j, k, :) = [rnd_cp innz abs_cp outnz abs_rspp rnd_rsp];
        end
    end
end

elapsed_time = toc;
fprintf('Finished, Time Taken %.4f Seconds\n', elapsed_time);

%% Clear Data
clear i j k Decompositions FcnValThresh MaxOuterItr NumItr Rank Rs Rc t_sz Tensor thresh T
clear elapsed_time abs_cp abs_rspp innz outnz rnd_cp rnd_rsp RSP_K SP_K K