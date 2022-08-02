function sm_ConvertKilosort3Neurosuite(rez)

% Converts KiloSort templates Klusta into klusters-compatible
% fet,res,clu,spk files.  Works on a single shank of a recording, assumes a
% 16bit .dat and an .xml file is present in "basepath" (home folder) and
% that they are named basename.dat and basename.xml.
%
% Inputs:
%   basepath -  directory path to the main recording folder with .dat and .xml
%               as well as shank folders made by makeProbeMapKlusta2.m (default is
%               current directory matlab is pointed to)
%   basename -  shared file name of .dat and .xml (default is last part of
%               current directory path, ie most immediate folder name)

% Brendon Watson 2016
% Edited by Peter Petersen 2017

savepath = rez.ops.savepath;
basepath = rez.ops.basepath;
basename = rez.ops.basename;
if ~exist('basepath','var')
    [~,basename] = fileparts(cd);
    basepath = cd;
end
if ~exist('rez','var');
    load(fullfile(basepath,'rez.mat'))
end
mkdir(fullfile(savepath,'OriginalClus'))


%%
Nchan =rez.ops.NchanTOT;
connected    = rez.connected;
%xcoords      = rez.xc;
%ycoords      = rez.yc;
% Nchan = rez.ops.Nchan;
% connected    = ones(Nchan, 1);
% xcoords      = ones(Nchan, 1);
% ycoords      = (1:Nchan)';

par = LoadXml(fullfile(basepath,[basename '.xml']));

totalch = par.nChannels;
sbefore = 24;%samples before/after for spike extraction
safter = 24;%... could read from SpkGroups in xml



% if isfield(par.SpkGrps,'nSamples')
%     if ~isempty(par.SpkGrps(1).nSamples);
%         if isfield(par.SpkGrps,'PeakSample')
%             if ~isempty(par.SpkGrps(1).PeakSample);
%                 sbefore = par.SpkGrps(1).PeakSample;
%                 safter = par.SpkGrps(1).nSamples - par.SpkGrps(1).PeakSample;
%             end
%         end
%     end
% end

if exist(rez.ops.fbinary,'file')
    datpath = rez.ops.fbinary;
else
    datpath = fullfile(basepath,[basename '.dat']);
end
%% spike extraction from dat
if ~exist('dat')
    dat             = memmapfile(datpath,'Format','int16');
    
end
% [spikeTimes, ii] = sort(spikeTimes);


spktimes = uint64(rez.ts);
clu = uint32(rez.clu);
clu = clu+1; %base 1

%%

%do a little QC

%only take units with >0.05Hz FR
fs = rez.ops.fs;

recLen = (length(dat.data)/totalch)/fs;


ts = double(spktimes)/fs ;
uclu = unique(clu);
[~,b] = ismember(clu,uclu);
n = histc(b,1:max(b));
kpu = (n/ts(end)) > .05 &  (n/ts(end)) < 60; % rate>.05 Hz and <100Hz

kp = ismember(clu,uclu(kpu)) & (ts+.002)< recLen ;

spktimes = spktimes(kp);
clu = clu(kp);




%% do homework for assigning templates to shanks
% [~,shank]=fileparts(basepath);
templates = rez.W;
% m = min(templates,[],2);%find the min value of each waveform on each channel
% [~,m] = min(m,[],1);%find which channel minimum is least
% m = squeeze(m);%which channel is minimum on each template
m = min(templates,[],2);%find the most deviated value of each waveform on each channel
m = squeeze(m);
[~,m] = min(m,[],2);%find which channel has most deviated value for each templnate


grouplookup = rez.ops.kcoords;%list of group/shank of each channel
templateshankassignments = grouplookup(m);%for the list of maximal channels, which group is each in
allgroups = unique(grouplookup);

%Grp 0 contain discared channels
allgroups(allgroups==0) = [];

for groupidx = 3:length(allgroups)
    
    %if isfield(par.SpkGrps(groupidx),'Channels')
    %if ~isempty(par.SpkGrps(groupidx).Channels)
    % for each group loop through, find all templates clus
    tgroup          = allgroups(groupidx);%shank number
    ttemplateidxs   = find(templateshankassignments==tgroup);%which templates/clusters are in that shank
    ttemplates      = squeeze(templates(ttemplateidxs,:,:));
    
    
    tidx            = ismember(clu,ttemplateidxs);%find spikes indices in this shank
    tclu            = clu(tidx);%extract template/cluster assignments of spikes on this shank
    tspktimes       = spktimes(tidx);
    
    
    %  tclu = tclu(maxSpk_ix);
    
    %  tspktimes = tspktimes(maxSpk_ix);
    
    
    
    gidx            = find(rez.ops.kcoords == tgroup);%find all channels in this group
    channellist     = [];
    
    for ch = 1:length(par.SpkGrps)
        if ismember(gidx(1),par.SpkGrps(ch).Channels+1)
            channellist = par.SpkGrps(ch).Channels+1;
            break
        end
    end
    if isempty(channellist)
        disp(['Cannot find spkgroup for group ' num2str(groupidx) ])
        continue
    end
    
    
    tsampsperwave   = (sbefore+safter);
    ngroupchans     = length(channellist);
    %  valsperwave     = tsampsperwave/4 * ngroupchans;
    valsperwave     = tsampsperwave * ngroupchans;
    %  wvforms_all     = zeros(length(tspktimes)*tsampsperwave/4*ngroupchans,1,'int16');
    
    wvforms_all     = zeros(length(tspktimes)*tsampsperwave*ngroupchans,1,'int16');
    
    
    wvranges        = zeros(length(tspktimes),ngroupchans);
    wvpowers        = zeros(1,length(tspktimes));
    
    for j=1:length(tspktimes)
        
        w = dat.data((double(tspktimes(j))-sbefore).*totalch+1:(double(tspktimes(j))+safter).*totalch);
        wvforms=reshape(w,totalch,[]);
        %select needed channels
        wvforms = wvforms(channellist,:);
        %         % detrend
        %         wvforms = floor(detrend(double(wvforms)));
        % median subtract
        wvforms = wvforms - repmat(median(wvforms')',1,sbefore+safter);
        %  wvforms = wvforms(:,1:4:end);%downsample to
        wvforms = wvforms(:);
        
        
        %some processing for fet file
        wvaswv = reshape(wvforms,tsampsperwave,ngroupchans);
        wvranges(j,:) = range(wvaswv);
        wvpowers(j) = sum(sum(wvaswv.^2));
        
        %  lastpoint = tsampsperwave/4*ngroupchans*(j-1);
        
        lastpoint = tsampsperwave*ngroupchans*(j-1);
        wvforms_all(lastpoint+1 : lastpoint+valsperwave) = wvforms;
        %     wvforms_all(j,:,:)=int16(floor(detrend(double(wvforms)')));
        if rem(j,100000) == 0
            disp([num2str(j) ' out of ' num2str(length(tspktimes)) ' done'])
        end
    end
    wvranges = wvranges';
    
    
    %%
    if ~isempty(wvforms_all)
        
        %spk save
        savepath_fma =  rez.ops.basepath;
        spkname = fullfile(savepath_fma, [basename '.spk.' num2str(tgroup)]);
        
        wvforms_all = reshape(wvforms_all,ngroupchans,tsampsperwave,[]);
        waveforms1 = wvforms_all(:,1:4:end,:); % downsample
        waveforms1 = waveforms1(:);
        
        
        
        
        fid=fopen(spkname,'w');
        fwrite(fid,waveforms1,'int16');
        fclose(fid);
        clear fid
        
        
        %% Spike features
        
        
        fetname = fullfile(savepath_fma, [basename '.fet.' num2str(tgroup)]);
        
        nComp = 3;
        fet = nan(size(wvforms_all,3),nComp*size(wvforms_all,1)+2);
        for i = 1:size(wvforms_all,1)
            e_d = squeeze(wvforms_all(i,:,:));
            [~, score] = pca(zscore(double(e_d')));
            
            fet(:,(i-1)*nComp+1:(i-1)*nComp+nComp) = score(:,1:nComp);
        end
        
        fet(:,end-1) = wvpowers;
        fet(:,end) = double(tspktimes);
        
        
        
        SaveFetIn(fetname,fet);
        
        
        %%
        cluname = fullfile(savepath_fma, [basename '.clu.' num2str(tgroup)]);
        tclu = cat(1,length(unique(tclu)),double(tclu));
        fid=fopen(cluname,'w');
        %     fprintf(fid,'%d\n',clu);
        fprintf(fid,'%.0f\n',tclu);
        fclose(fid);
        clear fid
        %%
        
        %res
        resname = fullfile(savepath_fma, [basename '.res.' num2str(tgroup)]);
        fid=fopen(resname,'w');
        fprintf(fid,'%.0f\n',tspktimes);
        fclose(fid);
        clear fid
        
        
        %%
        disp(['Shank ' num2str(tgroup) ' done'])
        
        
        
        
    end
    %copyfile(fullfile(savepath_fma, [basename,'.clu.*']),fullfile(savepath, 'OriginalClus'))
end
clear dat

end

function SaveFetIn(FileName, Fet, BufSize)

if nargin<3 | isempty(BufSize)
    BufSize = inf;
end

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
    formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');
fprintf(outputfile, '%d\n', nFeatures);

if isinf(BufSize)
    
    temp = [round(100* Fet(:,1:end-1)) round(Fet(:,end))];
    fprintf(outputfile,formatstring,temp');
else
    nBuf = floor(size(Fet,1)/BufSize)+1;
    
    for i=1:nBuf
        BufInd = [(i-1)*nBuf+1:min(i*nBuf,size(Fet,1))];
        temp = [round(100* Fet(BufInd,1:end-1)) round(Fet(BufInd,end))];
        fprintf(outputfile,formatstring,temp');
    end
end
fclose(outputfile);
end