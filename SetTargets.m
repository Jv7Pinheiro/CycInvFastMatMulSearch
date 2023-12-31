function target = SetTargets(Mat, thresh)
    % Set target to previous solution (for penalty term)
    target = Mat;
    % Cap large values
    target = min(target,1);
    target = max(target,-1);

    % Set small entries to zero
    mask = abs(target) > thresh;
    target = target .* mask;

    target = roundWithThreshold(target, thresh);
end