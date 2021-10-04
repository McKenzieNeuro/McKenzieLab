
function sm_Process_Intan_kilosort(fbasename,dirName,varargin)

% Process_Intan(fbasename,mergename)
% Process the Intan data folders starting with 'fbasename', rename each
% recording 'fbasename-01,2,3' in chronological order and lauch
% process_multi_start


cd(dirName)





% parse args
p = inputParser;
addParameter(p,'mergename',[],@isstr)
addParameter(p,'intan_version','Combined',@isstr)



parse(p,varargin{:})
mergename = p.Results.mergename;
intan_version = p.Results.intan_version;




if strcmpi(intan_version,'OpenEphys')
    basename = dl_BasenameFromBasepath(dirName);
else
    basename = bz_BasenameFromBasepath(dirName);
end

if isempty(mergename)
    [~,mergename,~] = fileparts(fbasename);
    
end


fprintf('Processing %s...\n',mergename);

%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
% convert file format
%%%%%%%%%%%%%%%%%%%%%

recList = sm_Process_IntanData_multi_start(fbasename,'intan_version',intan_version);

%maybe already ran clean up, check for files in right format
if isempty(recList)
    datf = getAllExtFiles(pwd,'dat',1);
    kp = cellfun(@any,regexp(datf,[mergename '-[0-9]+.dat']));
    recList =  datf(kp);
    for ii=1:length(recList)
        recList{ii} = recList{ii}(1:end-4);
    end
    
end

%if ~ isempty(recList)
%    UpdateXml_MergeName([mergename '-01'],mergename);
%end



xmlfil = getAllExtFiles(dirName,'xml',0);
%% concatenate files
datfil = getAllExtFiles(dirName,'dat',0);
datfil = datfil(cellfun(@any,regexp(datfil,[basename '-[0-9]'])));


digitalfils = cellfun(@any,regexp(datfil,'digitalin'));

digfils = datfil(digitalfils);

outfil_xml = [basename '.xml'];
outfil_dat = [basename '.dat'];
digitalfils_dat = [basename '_digitalin.dat'];

if strcmpi(intan_version,'combined')
    
    analogfils = [];
    auxfils = [];
    
    datafil = ~ (digitalfils);
    
    datafil = datfil(datafil);
    
    
else
    %handle analogue and auxiliary channels
    
    analogfils = cellfun(@any,regexp(datfil,'analogin'));
    auxfils = cellfun(@any,regexp(datfil,'auxiliary'));
    
    auxfils = datfil(auxfils);
    analogfils = datfil(analogfils);
    
    analogfils_dat = [basename '_analogin.dat'];
    auxfils_dat = [basename '_auxiliary.dat'];
    
    datafil = ~ (analogfils | auxfils | digitalfils);
    
    datafil = datfil(datafil);
    
end


%data file
if length(datafil)>1
    cmd = [];
    for i = 1:length(datafil)
        cmd = [cmd datafil{i} ' ' ];
    end
    cmd = ['type ' cmd ' > ' outfil_dat];
    status = system(cmd);
    
    if status ==0
        for i = 1:length(datafil)
            delete(datafil{i} )
            
            if i ==1
                
                movefile(xmlfil{i} , outfil_xml);
                
            else
                delete(xmlfil{i} )
                
            end
            
        end
    end
    
    
    
elseif length(datafil)==1
    movefile( datafil{1} , outfil_dat);
    
    if exist( xmlfil{1})
        movefile( xmlfil{1}, outfil_xml);
    end
end


%analogin file
if length(analogfils)>1
    cmd = [];
    for i = 1:length(analogfils)
        cmd = [cmd analogfils{i} ' ' ];
    end
    cmd = ['type ' cmd ' > ' analogfils_dat];
    status = system(cmd);
    
    if status ==0
        for i = 1:length(analogfils)
            delete( analogfils{i} );
        end
    end
    
    
elseif  length(analogfils)==1
    movefile( analogfils{1} , analogfils_dat);
    
end



%digitalin file
if length(digfils)>1
    cmd = [];
    for i = 1:length(digfils)
        cmd = [cmd digfils{i} ' ' ];
    end
    cmd = ['type ' cmd ' > ' digitalfils_dat];
    status = system(cmd);
    
    if status ==0
        for i = 1:length(digfils)
            delete( digfils{i} );
            
        end
    end
    
    
elseif  length(digfils)==1
    movefile(digfils{1} , digitalfils_dat);
    
end


if length(auxfils)>1
    cmd = [];
    for i = 1:length(auxfils)
        cmd = [cmd auxfils{i} ' ' ];
    end
    cmd = ['type ' cmd ' > ' auxfils_dat];
    status = system(cmd);
    
    if status == 0
        for i = 1:length(auxfils)
            delete(auxfils{i} );
            
        end
    end
    
elseif  length(auxfils)==1
    movefile( auxfils{1} , auxfils_dat);
    
end





end
