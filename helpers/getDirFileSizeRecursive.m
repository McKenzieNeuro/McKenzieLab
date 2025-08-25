function totalSize = getDirFileSizeRecursive(dirPath)
% getDirFileSizeRecursive Recursively calculates total size of all files
%
%   totalSize = getDirFileSizeRecursive(dirPath)
%
%   Input:
%       dirPath - root directory path (string or char)
%   Output:
%       totalSize - total size in bytes of all files under this directory

    totalSize = 0;

    % Get all items in the current directory
    entries = dir(dirPath);

    % Loop through all items
    for i = 1:length(entries)
        entry = entries(i);

        % Skip '.' and '..'
        if strcmp(entry.name, '.') || strcmp(entry.name, '..')
            continue;
        end

        fullPath = fullfile(dirPath, entry.name);

        if entry.isdir
            % Recurse into subdirectory
            totalSize = totalSize + getDirFileSizeRecursive(fullPath);
        else
            % Add file size
            totalSize = totalSize + entry.bytes;
        end
    end
end