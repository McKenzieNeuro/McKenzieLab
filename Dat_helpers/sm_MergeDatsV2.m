function sm_MergeDats(fnameIn, fnameOut, channels, nbChan)
% sm_MergeDats merges multiple temporary .dat files into a single binary.
%   fnameIn   = cell array of input filenames
%   fnameOut  = output filename
%   channels  = row index or indices to extract from each file (default = 1)
%   nbChan    = number of channels per input file (default = 1)

    % Handle default arguments
    if nargin < 3
        channels = 1;
    end
    if nargin < 4
        nbChan = 1;
    end

    sizeInBytes = 2;   % bytes per int16 sample
    chunk       = 1e5; % samples per chunk

    % Determine total number of elements in the first file
    fInfo      = dir(fnameIn{1});
    nElements  = fInfo.bytes / sizeInBytes;
    chunkSize  = chunk * nbChan;
    nbChunks   = floor(nElements / chunkSize);

    % Open output file
    fidO = fopen(fnameOut, 'w');
    nFiles = numel(fnameIn);
    fidI   = zeros(1, nFiles);
    for iFile = 1:nFiles
        fidI(iFile) = fopen(fnameIn{iFile}, 'r');
    end

    % Process main chunks
    for ii = 1:nbChunks
        dat = nan(nFiles, chunk);
        for jj = 1:nFiles
            tmp = fread(fidI(jj), chunkSize, 'int16');
            tmp = reshape(tmp, nbChan, []);
            dat(jj, :) = tmp(channels, :);
        end
        fwrite(fidO, dat(:), 'int16');
    end

    % Process remainder samples without large allocation
    remainderElems = mod(nElements, chunkSize);
    if remainderElems > 0
        for jj = 1:nFiles
            tmp = fread(fidI(jj), remainderElems, 'int16');
            tmp = reshape(tmp, nbChan, []);
            fwrite(fidO, tmp(channels, :), 'int16');
        end
    end

    % Cleanup file handles
    fclose(fidO);
    for iFile = 1:nFiles
        fclose(fidI(iFile));
    end
end