%% Set Up
MMT2 = matmul_tensor(2, 2, 2); 
% Rank 7 | 2 Decompositions
% (Rs = 4, Rc = 1), (Rs = 1, Rc = 2).

MMT3 = matmul_tensor(3, 3, 3);
% Rank 23 | 7 Decompositions
% (Rs = 20, Rc = 1), (Rs = 17, Rc = 2), (Rs = 14, Rc = 3),
% $(Rs = 11, Rc = 4)$, (Rs = 8, Rc = 5), $(Rs = 5, Rc = 6)$, 
% $(Rs = 2, Rc = 7)$
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
% (Rs = 19, Rc = 10), $(Rs = 16, Rc = 11)$, (Rs = 13, Rc = 12), 
% (Rs = 10, Rc = 13), (Rs = 7, Rc = 14), (Rs = 4, Rc = 15), 
% $(Rs = 1, Rc = 16)$
% Rank 48 | 15 Decompositions
% (Rs = 45, Rc = 1), (Rs = 42, Rc = 2), (Rs = 39, Rc = 3),
% (Rs = 36, Rc = 4), (Rs = 33, Rc = 5), (Rs = 30, Rc = 6),
% (Rs = 27, Rc = 7), (Rs = 24, Rc = 8), (Rs = 21, Rc = 9),
% (Rs = 18, Rc = 10), (Rs = 15, Rc = 11), (Rs = 12, Rc = 12),
% (Rs = 9, Rc = 13), (Rs = 6, Rc = 14), (Rs = 3, Rc = 15),

MMT5 = matmul_tensor(5, 5, 5);
% Rank 99 | 32 Decompositions
% (Rs = 96, Rc = 1), (Rs = 93, Rc = 2), (Rs = 90, Rc = 3), 
% (Rs = 87, Rc = 4), (Rs = 84, Rc = 5), (Rs = 81, Rc = 6), 
% (Rs = 78, Rc = 7), (Rs = 75, Rc = 8), (Rs = 72, Rc = 9), 
% (Rs = 69, Rc = 10), (Rs = 66, Rc = 11), (Rs = 63, Rc = 12), 
% (Rs = 60, Rc = 13), (Rs = 57, Rc = 14), (Rs = 54, Rc = 15), 
% (Rs = 51, Rc = 16), (Rs = 48, Rc = 17), (Rs = 45, Rc = 18), 
% (Rs = 42, Rc = 19), (Rs = 39, Rc = 20), (Rs = 36, Rc = 21), 
% (Rs = 33, Rc = 22), (Rs = 30, Rc = 23), (Rs = 27, Rc = 24), 
% (Rs = 24, Rc = 25), (Rs = 21, Rc = 26), (Rs = 18, Rc = 27), 
% (Rs = 15, Rc = 28), (Rs = 12, Rc = 29), (Rs = 9, Rc = 30), 
% (Rs = 6, Rc = 31), (Rs = 3, Rc = 32)
% Rank 98 | 32 Decompositions
% (Rs = 95, Rc = 1), (Rs = 92, Rc = 2), (Rs = 89, Rc = 3)
% (Rs = 86, Rc = 4), (Rs = 83, Rc = 5), (Rs = 80, Rc = 6)
% (Rs = 77, Rc = 7), (Rs = 74, Rc = 8), (Rs = 71, Rc = 9)
% (Rs = 68, Rc = 10), (Rs = 65, Rc = 11), (Rs = 62, Rc = 12)
% (Rs = 59, Rc = 13), (Rs = 56, Rc = 14), (Rs = 53, Rc = 15)
% (Rs = 50, Rc = 16), (Rs = 47, Rc = 17), (Rs = 44, Rc = 18)
% (Rs = 41, Rc = 19), (Rs = 38, Rc = 20), (Rs = 35, Rc = 21)
% (Rs = 32, Rc = 22), (Rs = 29, Rc = 23), (Rs = 26, Rc = 24)
% (Rs = 23, Rc = 25), (Rs = 20, Rc = 26), (Rs = 17, Rc = 27)
% (Rs = 14, Rc = 28), (Rs = 11, Rc = 29), (Rs = 8, Rc = 30)
%(Rs = 5, Rc = 31), (Rs = 2, Rc = 32)
% Rank 97 | 32 Decompositions
% (Rs = 94, Rc = 1), (Rs = 91, Rc = 2), (Rs = 88, Rc = 3), 
% (Rs = 85, Rc = 4), (Rs = 82, Rc = 5), (Rs = 79, Rc = 6), 
% (Rs = 76, Rc = 7), (Rs = 73, Rc = 8), (Rs = 70, Rc = 9), 
% (Rs = 67, Rc = 10), (Rs = 64, Rc = 11), (Rs = 61, Rc = 12), 
% (Rs = 58, Rc = 13), (Rs = 55, Rc = 14), (Rs = 52, Rc = 15), 
% (Rs = 49, Rc = 16), (Rs = 46, Rc = 17), (Rs = 43, Rc = 18), 
% (Rs = 40, Rc = 19), (Rs = 37, Rc = 20), (Rs = 34, Rc = 21), 
% (Rs = 31, Rc = 22), (Rs = 28, Rc = 23), (Rs = 25, Rc = 24), 
% (Rs = 22, Rc = 25), (Rs = 19, Rc = 26), (Rs = 16, Rc = 27), 
% (Rs = 13, Rc = 28), (Rs = 10, Rc = 29), (Rs = 7, Rc = 30), 
% (Rs = 4, Rc = 31), (Rs = 1, Rc = 32)

% Set which tensor to test
T = MMT5; % Decomposing Tensor

if isequal(T, MMT2)
    NumItr = 50;
    Tensor = 'MMT2';
    
    Rs = 1;
    Rc = 2;
    
    clear MMT3 MMT4;
elseif isequal(T, MMT3)
    NumItr = 10000;
    Tensor = 'MMT3';

    Rs = 8;
    Rc = 5;

    clear MMT2 MMT4;
elseif isequal(T, MMT4)
    NumItr = 10000;
    Tensor = 'MMT4';
    
    Rs = 21;
    Rc = 9;

    clear MMT2 MMT3;
elseif isequal(T, MMT5)
    NumItr = 15000;
    Tensor = 'MMT5';

    Rs = 42;
    Rc = 19;

    clear MMT2 MMT3 MMT4
end

% Thresholds used in CI_sparsify, roundWithThresholds, and Setting Targets
thresh = [0.01 0.05 0.1 0.2 0.3 0.4 0.5]';
t_sz = size(thresh, 1);

% This is how many times we will use a out_cell as in_cell
MaxOuterItr = 25;

%% Start Testing
Data = zeros(NumItr, t_sz, MaxOuterItr, 3);
Matrices = cell(size(Data(:, :, :, 1)));
fprintf('\nSearching %s solutions with Rs=%d, Rc=%d\n', Tensor, Rs, Rc);
fprintf('Settings:\n     %d Number of iterations\n     %d Outer Iterations\n     Thresholds:\n', NumItr, MaxOuterItr);
disp(thresh)


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
                [K, ~, out] = ci_cp_dgn(T, Rs, Rc, 'printitn', 0, 'maxiters', 150, 'lambda', 1e-6, 'tol', 1e-8);
            else
                [K, ~, out] = ci_cp_dgn(T, Rs, Rc, 'init', RSP_K, 'printitn', 0, 'maxiters', 150, 'lambda', 1e-6, 'tol', 1e-8);
            end
            [rnd_cp, ~, abs_cp] = getErrors(K);
            [SP_K] = CI_sparsify(K, thresh(j));

            % Rounding with Threshold
            RSP_K = cellfun(@(x) roundWithThreshold(x, thresh(j)), SP_K, 'UniformOutput', false); 
             
            Data(i, j, k, :) = [out.FcnVal rnd_cp abs_cp];
            Matrices{i, j, k} = K;
        end
    end
end

elapsed_time = toc;
fprintf('Finished, Time Taken %.4f Seconds\n', elapsed_time);
clear T Tensor
save("Data_99_42_19");
%% Clear Data
clear i j k Decompositions FcnValThresh MaxOuterItr NumItr Rank Rs Rc t_sz Tensor thresh T
clear abs_cp abs_rspp innz outnz rnd_cp rnd_rsp RSP_K SP_K K