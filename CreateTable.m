function T = CreateTable(Data)
    Rows = length(Data);
    Columns = 1 +numel(fieldnames(Data{1}.structData));
    
    T = table(1, 2, 3, 4, 5, 6, 'VariableNames', {'rng', 'ExitFlag', 'FcnVal', 'RelErr', 'Fit', 'ExitDiff'});

    for i = 1:Rows
        newRow = struct2table(Data{i}.structData);
        newColumn = array2table(Data{i}.numericData, "VariableNames", {'rng'});
        combinedRow = [newColumn, newRow];

        T = [T; combinedRow];
    end

    T(1, :) = [];
end