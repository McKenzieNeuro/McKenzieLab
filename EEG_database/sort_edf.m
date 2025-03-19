fils=  getAllExtFiles('N:\Research-Studies\Study 21-044 HRDI BM','edf',1);
%%
out = 'N:\Research-Studies\Study 21-044 HRDI BM\deidentified\edf';

for i = 1:length(fils)
    
    outdir = [out filesep fils{i}(42:end)];
    
    if ~exist(fileparts(outdir))
        sm_makeNewDir(fileparts(outdir))
    end
    
    copyfile(fils{i},outdir)
    
end

%%
fils=  getAllExtFiles('N:\Research-Studies\Study 21-044 HRDI BM\deidentified','edf',1);
%%
ii=1;
for i = 2:length(fils)
    
    try
        [header] = edfread2(fils{i});
        % get old day
        Y = str2num(['20' header.startdate(7:8)]);
        M = str2num(header.startdate(4:5));
        D =  str2num(header.startdate(1:2));
        
        sl = regexp(fils{i},filesep);
        subjID = fils{i}(sl(5)+1:sl(6)-1);
        subjID = strrep(subjID,' ','');
        subjIDo = ['X F X ' subjID];
        t = datetime(Y,M,D) - caldays(7);
        yr = num2str(t.Year) ;
        yr = yr(3:4);
        
        
        
        
        
        % define new out dir
        outdir = 'N:\Research-Studies\Study 21-044 HRDI BM\deidentified\edf1\';
        subdir = fils{i}(59:end);
        
        % deal with subdir
        [subdir,b,c]  =fileparts(subdir);
        
        % if it is scored as patient
        sd = regexp(subdir,'Patient');
        
        
        
        sl = regexp(subdir,filesep);
        
        for j = 1:length(sd)
            
            
            idx1 = sl(find(sl>sd(j),1,'first'));
            if j<length(sd)
               
                idx2 = sd(j+1)-1;
            else
                
                idx2 = length(subdir);
            end
            
            if j==1
                newsubdir = subdir(1:sd(j)-1);
                newsubdir = [newsubdir 'Patient_' datestr(t)];
            else
                newsubdir = [newsubdir 'Patient_' datestr(t)];
                
            end
            
            newsubdir = [newsubdir subdir(idx1:idx2)];
        end
        
        newFname = [subjID '_' datestr(t) '_'  strrep(header.starttime,'.','-')];
        outfil = [outdir newsubdir filesep newFname '.edf'];
         if ~exist(outfil)
            if ~exist(fileparts(outfil))
                sm_makeNewDir(fileparts(outfil))
            end

         end
        cmd = ['C:\Users\samckenzie\PycharmProjects\deid_edf\main.py ' '''' fils{i} '''' ' '''  outfil ''''];
       
          if ~exist(outfil)
            if ~exist(fileparts(outfil))
                sm_makeNewDir(fileparts(outfil))
            end
             pyrunfile(cmd)
          end
     
    catch
        
        badfil{ii} = fils{i}
        save('N:\Research-Studies\Study 21-044 HRDI BM\deidentified\edf1\doOver.mat','badfil')
        ii = ii+1;
    end
end

%%

%
%      fid = fopen(fils{i}, 'r+');
%     fprintf(fid, '%-8s', '0');   % Version must be 0
% fprintf(fid, '%-80s', 'DEIDENTIFIED'); % Remove patient info
% fprintf(fid,'%-80s', 'DEIDENTIFIED'); %Remove recording info
% fprintf(fid, '%02i.%02i.%02i', 1,1,1); % Set date as 01.01.01
%
% fclose(fid);
