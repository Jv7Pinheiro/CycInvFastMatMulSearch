% Group Search attempts to search for solution of cyclic invariant 
% matrices per matmul_tensor

% The rank of each tensor is fixed, and to find what Rs and Rc are we
% use the function Rank = Rs + 3Rc

MMT2 = matmul_tensor(2, 2, 2);
MMT3 = matmul_tensor(3, 3, 3);
MMT4 = matmul_tensor(4, 4, 4);


%% MatMul_Tensor 4x4x4
savepath = "C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\GroupSearches\DataMMT2.mat";

% How many times the above tensor will decompose
NumItr = 500;

% Tensor's rank
Rank = 7;

% Searches are 
% (Rs = 4, Rc = 1), (Rs = 1, Rc = 2).

DataMMT2 = struct([]);
currData = struct([]);

% Iterate through Rc as it increases by one
fprintf('\nSearching MMT2 solutions\n')
for Rc = 1:2
    Rs = Rank - 3*Rc;
    
    SucItrNum = 1;

    for i = 1:NumItr
        rng(i);
        [K, P0, output] = ci_cp_dgn(MMT2, Rs, Rc, printitn = 0);
        if (output.FcnVal < 1.0e-04 && output.ExitDiff < 1)
            iterationData.numericData = i;
            iterationData.cellData = K;
            iterationData.structData = output;
    
            currData{SucItrNum} = iterationData;
            SucItrNum = SucItrNum + 1;
        end
    end
    DataMMT2{Rc} = currData';
    fprintf('Done with Rs=%d, Rc=%d trials\n', Rs, Rc);
end

DataMMT2 = DataMMT2';
fprintf('Done with MMT2 trials\n\n');
save(savepath, "DataMMT2");
%% MatMul_Tensor 9x9x9
savepath = "C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\GroupSearches\DataMMT3HalfLambda.mat";

% How many times tensor MMT3 will decompose
NumItr = 2000;

% Tensor's rank
Rank = 23;
% Searches are
% (Rs = 20, Rc = 1), (Rs = 17, Rc = 2), (Rs = 14, Rc = 3),
% (Rs = 11, Rc = 4), (Rs = 8, Rc = 5), (Rs = 5, Rc = 6), 
% (Rs = 2, Rc = 7)

DataMMT3 = struct([]);
currData = struct([]);

% Iterate through Rc as it increases by one
fprintf('\nSearching MMT3 solutions\n')
for Rc = 1:7
    Rs = Rank - 3*Rc;
    
    SucItrNum = 1;

    for i = 1:NumItr
        rng(i);
        [K, P0, output] = ci_cp_dgn(MMT3, Rs, Rc, printitn = 0, maxiters = 100, lambda = 1e-2);
        if (output.FcnVal < 1.0e-01 && output.ExitDiff < 1)
            iterationData.numericData = i;
            iterationData.cellData = K;
            iterationData.structData = output;
    
            currData{SucItrNum} = iterationData;
            SucItrNum = SucItrNum + 1;
        end
    end
    DataMMT3{Rc} = currData';
    fprintf('Done with Rs=%d, Rc=%d trials\n', Rs, Rc);
end

DataMMT3 = DataMMT3';
fprintf('Done with MMT3 trials\n\n');
save(savepath, "DataMMT3");
%% MatMul_Tensor 16x16x16
savepath = "C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\GroupSearches\DataMMT4.mat";

% How many times the above tensor will decompose
NumItr = 3000;

%Tensor's Rank
Rank = 49;
% Searches are
% (Rs = 46, Rc = 1), (Rs = 43, Rc = 2), (Rs = 40, Rc = 3),
% (Rs = 37, Rc = 4), (Rs = 34, Rc = 5), (Rs = 31, Rc = 6), 
% (Rs = 28, Rc = 7), (Rs = 25, Rc = 8), (Rs = 22, Rc = 9),
% (Rs = 19, Rc = 10), (Rs = 16, Rc = 11), (Rs = 13, Rc = 12), 
% (Rs = 10, Rc = 13), (Rs = 7, Rc = 14), (Rs = 4, Rc = 15), 
% (Rs = 1, Rc = 16)

DataMMT4 = struct([]);
currData = struct([]);

fprintf('\nSearching MMT4 solutions\n')
for Rc = 1:16
    Rs = Rank - 3*Rc;

    SucItrNum = 1;

    for i = 1:NumItr
        rng(i);
        [K, P0, output] = ci_cp_dgn(MMT4, Rs, Rc, printitn = 0, maxiters = 200, lambda = 1e-8);
        if (output.ExitDiff < 1)
            iterationData.numericData = i;
            iterationData.cellData = K;
            iterationData.structData = output;
    
            currData{SucItrNum} = iterationData;
            SucItrNum = SucItrNum + 1;
        end
    end
    DataMMT4{Rc} = currData';
    fprintf('Done with Rs=%d, Rc=%d trials\n', Rs, Rc);
end

DataMMT4 = DataMMT4';
fprintf('Done with MMT4 trials\n\n');
save(savepath, "DataMMT4");