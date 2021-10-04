
function Process_Intan_kilosort(fbasename,dirName,varargin)

% Process_Intan(fbasename,mergename)
% Process the Intan data folders starting with 'fbasename', rename each
% recording 'fbasename-01,2,3' in chronological order and lauch
% process_multi_start


cd(dirName)
basename = bz_BasenameFromBasepath(dirName);
if isempty(varargin)
    [~,mergename,~] = fileparts(pwd);
else
    mergename = varargin{1};
end
fprintf('Processing %s...\n',mergename);

%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
% convert file format
%%%%%%%%%%%%%%%%%%%%%

recList = sm_Process_IntanData_multi_start(fbasename,mergename);

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
%%concatenate files
datfil = getAllExtFiles(dirName,'dat',0);

datfil = datfil(cellfun(@any,regexp(datfil,[basename '-[0-9]'])));
digitalfils = cellfun(@any,regexp(datfil,'digitalin'));
analogfils = cellfun(@any,regexp(datfil,'analogin'));
auxfils = cellfun(@any,regexp(datfil,'auxiliary'));
datafil = ~ (analogfils | auxfils | digitalfils);

datafil = datfil(datafil);
auxfils = datfil(auxfils);
digfils = datfil(digitalfils);
analogfils = datfil(analogfils);

outfil_xml = [basename '.xml'];
outfil_dat = [basename '.dat'];
analogfils_dat = [basename '_analogin.dat'];
auxfils_dat = [basename '_auxiliary.dat'];
digitalfils_dat = [basename '_digitalin.dat'];
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
  
    
    movefile( xmlfil{1}, outfil_xml);
   
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
