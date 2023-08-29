%% Testing Certain Settings

TCells3114 = cell(4, 3, 3, 3);
TOutputs3114 = cell(4, 3, 3, 3);
% We want to test, Sparcify Threshold, Set Targets Threshold, and Beta
TestsNum = [0.1, 0.05, 0.01]';
Tsz = numel(TestsNum);

TestingCells = {A_3114_1552 A_3114_1698 A_3114_567 A_3114_939}';
Csz = size(TestingCells, 1);

for i = 1:Csz % Iterate Through Desired Cells
    for j1 = 1:Tsz % Iterate Through Sparcify Threshold
        for j2 = 1:Tsz % Iterate Through Set Targets Threshold
            for j3 = 1:Tsz % Iterate Through Betas
                fprintf('i=%d, j1=%d, j2=%d, j3=%d\n', i, j1, j2, j3);
                [TCells3114{i, j1, j2, j3}, ~, TOutputs3114{i, j1, j2, j3}] = ci_cp_dgn(MMT3, 11, 4, 'init', CI_sparsify(TestingCells{i}, TestsNum(j1)), 'thresh', TestsNum(j2), 'beta', TestsNum(j3), 'printitn', 0);
            end
        end
    end
end

%% Interpreting the Data

% Convert to a cell array (if structs have varying fields)
CellData = TCells3114;
OutputData = TOutputs3114;
% Now you can access cellData{i, j, k, l} to get the individual struct

% Accessing values in the cell array (assuming the field exists)
RelErrCell = cellfun(@(x) x.RelErr, OutputData);
FcnValCell = cellfun(@(x) x.FcnVal, OutputData);
AbsErrCell = cellfun(@(x) x.AbsErr, OutputData);
