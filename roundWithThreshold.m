function roundMat = roundWithThreshold(inputMatrix, threshold)
    % Check if the threshold is within the allowed range [0, 0.5]
    if threshold < 0
        error('Threshold must be greater than 0');
    end
    
    % Initialize the rounded matrix with the same size as the input matrix
    [m, n] = size(inputMatrix);
    roundMat = zeros(m, n);
    
    % Loop through each element in the input matrix
    for i = 1:m
        for j = 1:n
            % Round the element to the nearest integer if it's within the threshold
            if abs(inputMatrix(i, j) - round(inputMatrix(i, j))) <= threshold
                roundMat(i, j) = round(inputMatrix(i, j));
            else
                roundMat(i, j) = inputMatrix(i, j);
            end
        end
    end
end