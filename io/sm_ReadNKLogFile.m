
function hdr = sm_ReadNKLogFile(LogFile)

% Open file
fid = fopen(LogFile, 'rb');
if (fid ~= -1)
   
    
    % Get file signature
    device = fread(fid, [1 16], '*char');
    if (get_header_version(device) == 0)
        error(['LOG file has unknown signature: ' device]);
    end
    % Get log blocks
    fseek(fid, 145, 'bof');
    n_logblocks = fread(fid, 1, 'uint8');
    % Initializations
    total_logs = 0;
    
    % Loop on log blocks
    for i = 1:n_logblocks
        % Read number of logs in this block
        fseek(fid, 146 + ((i-1) * 20) , 'bof');
        logblock_address = fread(fid, 1, 'uint32');
        fseek(fid, logblock_address + 18, 'bof');
        n_logs = fread(fid, 1, 'uint8');
        % Initialization
        fseek(fid, logblock_address + 20, 'bof');
        hdr.logs(i).label = cell(1, n_logs);
        hdr.logs(i).time  = zeros(1, n_logs);
        % Read all the events
        for j = 1:n_logs
            hdr.logs(i).label{j} = str_clean(fread(fid, [1 20], '*char'));
            timeH = str2double(fread(fid, [1 2], '*char'));
            timeM = str2double(fread(fid, [1 2], '*char'));
            timeS = str2double(fread(fid, [1 2], '*char'));
            hdr.logs(i).time(j) = 60*60*timeH + 60*timeM + timeS;
            hdr.logs(i).label2{j} = fread(fid, [1 19], '*char');
            % Compute time stamp
            timeH = str2double(hdr.logs(i).label2{j}(8:9));
            timeM = str2double(hdr.logs(i).label2{j}(10:11));
            timeS = str2double(hdr.logs(i).label2{j}(12:13));
            hdr.logs(i).timestamp(j) = 60*60*timeH + 60*timeM + timeS;
        end
        
        % Read sub-events
        try
            % Read number of sub-logs
            fseek(fid, 146 + (((i-1) + 22) * 20) , 'bof');
            sublogblock_address = fread(fid, 1, 'uint32');
            fseek(fid, sublogblock_address + 18, 'bof');
            n_sublogs = fread(fid, 1, 'uint8');
            % Read sub-logs
            if (n_sublogs == n_logs)
                fseek(fid, sublogblock_address + 20, 'bof');
                for j = 1:n_logs
                    hdr.logs(i).sublog{j} = fread(fid, [1 45], '*char');
                    hdr.logs(i).time(j) = hdr.logs(i).time(j) + str2double(['0.' hdr.logs(i).sublog{j}(25:30)]);
                end
            end
        catch
            disp('NK> Could not read sub-events.');
        end
        total_logs = total_logs + n_logs;
    end
    % Close file
    fclose(fid);
else
    hdr.logs = [];
end
end


%% ===== CHECK DEVICE =====
function ver = get_header_version(str)
    % Older NK systems
    if ismember(str, {...
        'EEG-1100A V01.00', ...
        'EEG-1100B V01.00', ...
        'EEG-1100C V01.00', ...
        'QI-403A   V01.00', ...
        'QI-403A   V02.00', ...
        'EEG-2100  V01.00', ...
        'EEG-2100  V02.00', ...
        'DAE-2100D V01.30', ...
        'DAE-2100D V02.00', ...
        'EEG-1100A V02.00', ...
        'EEG-1100B V02.00', ...
        'EEG-1100C V02.00'})
        ver = 1;
        
    % Newer NK systems (>= 2015)
    elseif ismember(str, {...
        'EEG-1200A V01.00'})
        ver = 2;
        
    else
        ver = 0;
    end
end

%% ===== CLEAN STRINGS =====
function s = str_clean(s)
    % Stop string at first termination
    iNull = find(s == 0, 1);
    if ~isempty(iNull)
        s(iNull:end) = [];
    end
    % Remove weird characters
    s(~ismember(s, '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-,:;.*+=?!<>''"`&%$()[]{}/\_@ áÁàÀâÂäÄãÃåÅæÆçÇéÉèÈêÊëËíÍìÌîÎïÏñÑóÓòÒôÔöÖõÕøØßúÚùÙûÛüÜ')) = [];
    % Remove useless spaces
    s = strtrim(s);
end

