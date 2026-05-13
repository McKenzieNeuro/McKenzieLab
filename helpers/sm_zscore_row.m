function zstack = sm_zscore_row(data)


% Get size
[rows, cols, slices] = size(data);

% Reshape tensor to 2D: combine slices into one larger 2D matrix
data_2d = reshape(data, rows, cols*slices);

% Z-score each row across columns
z_2d = zscore(data_2d, [], 2);

% Reshape back to original 3D shape
zstack = reshape(z_2d, rows, cols, slices);
end