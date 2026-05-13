function f1 = computeF1Score(y_true, y_pred)
    % Ensure inputs are column vectors
    y_true = y_true(:);
    y_pred = y_pred(:);

    % Calculate True Positives, False Positives, False Negatives
    TP = sum((y_true == 1) & (y_pred == 1));
    FP = sum((y_true == 0) & (y_pred == 1));
    FN = sum((y_true == 1) & (y_pred == 0));

    % Avoid division by zero
    if TP + FP == 0
        precision = 0;
    else
        precision = TP / (TP + FP);
    end

    if TP + FN == 0
        recall = 0;
    else
        recall = TP / (TP + FN);
    end

    % Compute F1 Score
    if precision + recall == 0
        f1 = 0;
    else
        f1 = 2 * (precision * recall) / (precision + recall);
    end
end