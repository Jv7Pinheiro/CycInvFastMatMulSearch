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
    Rank = 7;
    Decompositions = 2;
    NumItr = 10;
    Tensor = 'MMT2';

    FcnValThresh = 100;

    clear MMT3 MMT4;
elseif isequal(T, MMT3)
    Rank = 19;
    Decompositions = 6;
    NumItr = 5000;
    Tensor = 'MMT3';

    FcnValThresh = 1;

    clear MMT2 MMT4;
elseif isequal(T, MMT4)
    Rank = 49;
    Decompositions = 16;
    NumItr = 10000;
    Tensor = 'MMT4';

    FcnValThresh = 100000;

    clear MMT2 MMT3;
end

% Thresholds used in CI_sparsify, roundWithThresholds, and Setting Targets
thresh = [0.1 0.2 0.3 0.4 0.5]';
t_sz = size(thresh, 1);

% This vector is to measure how many succesful iterations there were per
% decomposition type
SucItrVec = zeros(Decompositions, 1);

% This is how many times we will use a out_cell as in_cell
MaxOuterItr = 7;
%% Testing
Prime_Vector = struct([]);
SP_Data = struct([]);
CP_Data = struct([]);


fprintf('\nSearching %s solutions with Rank %d\n', Tensor, Rank);

for Rc = 1:Decompositions
    Rs = Rank - 3*Rc;
    
    Curr_Prime_Vector = struct([]);
    Curr_SP_Data = struct([]);
    Curr_CP_Data = struct([]);

    SucItrNum = 1;
    
    tic;
    for i = 1:NumItr
        % Set seed to be the current attempt
        rng(i);
        
        % Run initial test with randomized values from seed
        [K_Prime, ~, output] = ci_cp_dgn(T, Rs, Rc, 'printitn', 0, 'maxiters', 150, 'lambda', 1e-6, 'tol', 1e-8);
        % Only save data if outputs are below thresholds
        if (output.FcnVal<FcnValThresh)
            [rnd, rel, abs] = getErrors(K_Prime);
            errors.rnd = rnd;
            errors.rel = rel;
            errors.abs = abs;
            

            iterationData1.num_val = i; % Save rng seed
            iterationData1.out_cell = K_Prime; % Save CP result
            iterationData1.out_struct = output; 
            iterationData1.errors = errors;

            
            Curr_Prime_Vector{SucItrNum} = iterationData1;
            for j = 1:t_sz
                for k = 1:MaxOuterItr
                    if (k == 1)
                        % Schur Sparsify
                        SP_K = CI_sparsify(K_Prime, thresh(j));
        
                        % Eigenvalue Sparsify
                        % P = CI_sparsify(K, thresholds(j), 'eig');
                    else
                        P = Curr_CP_Data{SucItrNum, j, k - 1}.out_cell;
                        % Schur Sparsify
                        SP_K = CI_sparsify(P, thresh(j));
        
                        % Eigenvalue Sparsify
                        % P = CI_sparsify(K, thresholds(j), 'eig');
                    end

                    % Rounding with Threshold
                    RSP_K = cellfun(@(x) roundWithThreshold(x, thresh(j)), SP_K, 'UniformOutput', false); 
                    

                    [rnd, rel, abs] = getErrors(RSP_K);
                    rsp_errors.rnd = rnd;
                    rsp_errors.rel = rel;
                    rsp_errors.abs = abs;
                    

                    iterationData2.RSP_K = RSP_K;
                    iterationData2.sp_thresh = thresh(j);
                    iterationData2.rsp_errors = rsp_errors;
                    
                    Curr_SP_Data{SucItrNum, j, k} = iterationData2;

                    
                    [K, P0, cp_output] = ci_cp_dgn(T, Rs, Rc, 'init', RSP_K, 'printitn', 0, 'maxiters', 150, 'lambda', 1e-6, 'tol', 1e-8);
                    
                    [rnd, rel, abs] = getErrors(K);
                    cp_errors.rnd = rnd;
                    cp_errors.rel = rel;
                    cp_errors.abs = abs;
                    
                    
                    iterationData3.in_cell = P0;
                    iterationData3.out_cell = K;
                    iterationData3.cp_output = cp_output;
                    iterationData3.sp_thresh = thresh(j);  
                    iterationData3.cp_errors = cp_errors;
   
                    Curr_CP_Data{SucItrNum, j, k} = iterationData3;
                end
            end
            % If there was a succesful iterations, increase successful
            % numer of iterations and try next rng value
            SucItrNum = SucItrNum + 1;
        end
    end
    SucItrVec(Rc) = SucItrNum - 1;
    SP_Data{Rc} = Curr_SP_Data;
    CP_Data{Rc} = Curr_CP_Data;
    Prime_Vector{Rc} = Curr_Prime_Vector';
    
    elapsed_time = toc;
    fprintf('Finished with Rs=%d, Rc=%d trials - Successful Iterations %d/%d - Time Taken %.4f Seconds\n', Rs, Rc, SucItrNum-1, NumItr, elapsed_time);
end


SP_Data = SP_Data';
CP_Data = CP_Data';
Prime_Vector = Prime_Vector';

fprintf('Done\n\n')

%% Clear Variables
clear i j k;
clear P0 K SP_K RSP_K output cp_output;
clear thresh t_sz s_thresh FcnValThresh round_thresh;
clear abs rel rnd;
clear Tensor;
clear NumItr MaxOuterItr;
clear elapsed_time;
clear Rank Decompositions;
clear Curr_SP_Data Curr_CP_Data Curr_Prime_Vector;
clear Rs Rc;
clear rsp_errors cp_errors errors;
clear SucItrNum;
clear iterationData1 iterationData2 iterationData3;
clear P K_Prime
clear SucItrVec MMT3
%% Save Outputs
save("General_Tests.mat", 'SP_Data', "CP_Data", "Prime_Vector")