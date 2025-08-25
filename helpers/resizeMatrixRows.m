function B = resizeMatrixRows(A, targetRows)
    % Resize a matrix A to have exactly targetRows rows
    % by duplicating rows while preserving their original order

    [nRows, nCols] = size(A);

    if nRows == targetRows
        B = A;  % No resizing needed
    elseif nRows > targetRows
        error('The input matrix has more rows than the target. This function only duplicates rows.');
    else
        % Compute how many times each row should be duplicated
        % Create an index list of rows to sample
        indices = linspace(1, nRows, targetRows);
        indices = round(indices); % Round to nearest row index

        % Ensure we don't exceed original bounds
        indices = min(max(indices, 1), nRows);

        B = A(indices, :);
    end
end