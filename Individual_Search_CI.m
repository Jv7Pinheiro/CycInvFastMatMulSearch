%% Individual Search Settings

% For each MatMulTensor we want to save the solutions,
% that happens when diff < 1 or when f < 1.0e-01 (so save them too)
% and the rng that was used to find said solutions.

MMT2 = matmul_tensor(2, 2, 2);
MMT3 = matmul_tensor(3, 3, 3);

NumItr = 1000;
% "C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data"

%% n = 2, Rs = 1, Rc = 2
AllCombinedData212 = struct([]);
SucItrNum = 1;

savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData212';
for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output212] = ci_cp_dgn(MMT2, 1, 2);
    if (output212.FcnVal < 1.0e-01 && output212.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output212;
        
        AllCombinedData212{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData212");
%% n = 2. Rs = 4, Rc = 1
AllCombinedData241 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData241';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    [K, P0, output241] = ci_cp_dgn(MMT2, 4, 1);
    if (output241.FcnVal < 1.0e-01 && output241.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output241;
        
        AllCombinedData241{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData241");
%% n = 3, Rs = 2, Rc = 7
AllCombinedData327 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData327';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output327] = ci_cp_dgn(MMT3, 2, 7);
    if (output327.FcnVal < 1.0e-01 && output327.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output327;
        
        AllCombinedData327{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData327");
%% n = 3, Rs = 11, Rc = 4
AllCombinedData3114 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData3114';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output3114] = ci_cp_dgn(MMT3,11, 4);
    if (output3114.FcnVal < 1.0e-01 && output3114.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output327;
        
        AllCombinedData3114{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData3114");
%% n = 3, Rs = 5, Rc = 6
AllCombinedData356 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData356';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output356] = ci_cp_dgn(MMT3, 5, 6);
    if (output356.FcnVal < 1.0e-01 && output356.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output356;
        
        AllCombinedData356{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData356");
%% n = 3, Rs = 7, Rc = 5
AllCombinedData375 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData375';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output375] = ci_cp_dgn(MMT3, 7, 5);
    if (output375.FcnVal < 1.0e-01 && output375.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output375;
        
        AllCombinedData375{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData375");
%% n = 3, Rs = 14, Rc = 3
AllCombinedData3143 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData3143';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output3143] = ci_cp_dgn(MMT3, 14, 3);
    if (output3143.FcnVal < 1.0e-01 && output3143.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output3143;
        
        AllCombinedData3143{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData3143");
%% n = 3, Rs = 17, Rc = 2
AllCombinedData3172 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData3172';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output3172] = ci_cp_dgn(MMT3, 17, 2);
    if (output3172.FcnVal < 1.0e-01 && output3172.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output3172;
        
        AllCombinedData3172{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData3172");
%% n = 3, Rs = 20, Rc = 1
AllCombinedData3201 = struct([]);
SucItrNum = 1;
savepath = 'C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\Projects\Strassen_Gradient&Jacobian\CI_Data\AllCombinedData3201';

for i = 1:NumItr
    fprintf("\n\nrng(%d): ", i);
    rng(i);
    [K, P0, output3201] = ci_cp_dgn(MMT3, 20, 1);
    if (output3201.FcnVal < 1.0e-01 && output3201.ExitDiff < 1)
        % Save the iteration (which is the same as the rng), and the outputs        
        iterationData.numericData = i;
        iterationData.cellData = K;
        iterationData.structData = output3201;
        
        AllCombinedData3201{SucItrNum} = iterationData;
        SucItrNum = SucItrNum + 1;
    end
end

save(savepath,"AllCombinedData3201");
