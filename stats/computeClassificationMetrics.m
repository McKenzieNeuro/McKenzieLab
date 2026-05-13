function [f1, sensitivity, specificity] = computeClassificationMetrics(y_true, y_pred)
    % Ensure inputs are column vectors
    y_true = y_true(:);
    y_pred = y_pred(:);

    % Calculate confusion matrix components
    TP = sum((y_true == 1) & (y_pred == 1));
    TN = sum((y_true == 0) & (y_pred == 0));
    FP = sum((y_true == 0) & (y_pred == 1));
    FN = sum((y_true == 1) & (y_pred == 0));

    % Calculate precision
    if TP + FP == 0
        precision = 0;
    else
        precision = TP / (TP + FP);
    end

    % Sensitivity (Recall)
    if TP + FN == 0
        sensitivity = 0;
    else
        sensitivity = TP / (TP + FN);
    end

    % Specificity (True Negative Rate)
    if TN + FP == 0
        specificity = 0;
    else
        specificity = TN / (TN + FP);
    end

    % F1 Score
    if precision + sensitivity == 0
        f1 = 0;
    else
        f1 = 2 * (precision * sensitivity) / (precision + sensitivity);
    end
end