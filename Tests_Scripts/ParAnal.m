fprintf("Number of random starts: %d", size(Data, 1));

FcnValues = Data(:, :, :, 1);
RndErrors = Data(:, :, :, 2);
AbsErrors = Data(:, :, :, 3);

FcnValues_Results = zeros(size(FcnValues, 3), size(FcnValues, 2), 2);
RndErrors_Results = zeros(size(RndErrors, 3), size(RndErrors, 2), 2);
AbsErrors_Results = zeros(size(AbsErrors, 3), size(AbsErrors, 2), 2);

for i = 1:size(FcnValues, 3)
    [FcnValues_Results(i, :, 1), FcnValues_Results(i, :, 2)] = min(FcnValues(:, :, i));
end
fprintf("\n\nFunction Values:\n")
disp(FcnValues_Results(:, :, 1))

[a, b] = min(FcnValues_Results(:, :, 1));
[minFV, d] = min(a);
idx = FcnValues_Results(b(d), d, 2);
fprintf("Minimum Function Value at (%d, %d, %d) is:\n", idx, d, b(d));
disp(minFV)



for i = 1:size(AbsErrors, 3)
    [RndErrors_Results(i, :, 1), RndErrors_Results(i, :, 2)] = min(RndErrors(:, :, i));
end
%fprintf("\n\nRound Errors:\n")
%disp(RndErrors_Results(:, :, 1))

[a, b] = min(RndErrors_Results(:, :, 1));
[minRND, d] = min(a);
idx = RndErrors_Results(b(d), d, 2);
fprintf("Minimum Round Error at (%d, %d, %d) is:\n", idx, d, b(d));
disp(minRND)




for i = 1:size(FcnValues, 3)
    [AbsErrors_Results(i, :, 1), AbsErrors_Results(i, :, 2)] = min(AbsErrors(:, :, i));
end
%fprintf("\n\nAbsolute Errors:\n")
%disp(AbsErrors_Results(:, :, 1))

[a, b] = min(AbsErrors_Results(:, :, 1));
[minABS, d] = min(a);
idx = AbsErrors_Results(b(d), d, 2);
fprintf("Minimum Absolute Error at (%d, %d, %d) is:\n", idx, d, b(d));
disp(minABS)


clear i
clear a b d idx minFV minRND minABS

%%

for i = 1:size(RndErrors, 3)
    I = find(RndErrors(:, :, i) == 0);
    if size(I, 1) ~= 0
        fprintf("Round\n")
        i
        disp(I)
    end
end

for i = 1:size(AbsErrors, 3)
    I = find(AbsErrors(:, :, i) == 0);
    if size(I, 1) ~= 0
        fprintf("Absolute\n")
        i
        disp(I)
    end
end

clear i I
