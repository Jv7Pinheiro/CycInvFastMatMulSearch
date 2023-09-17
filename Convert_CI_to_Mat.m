function T = Convert_CI_to_Mat(D)
    A = [D{1} D{2} D{3} D{4}];
    B = [D{1} D{4} D{2} D{3}];
    C = [D{1} D{3} D{4} D{2}];
    T = {A B C};
end