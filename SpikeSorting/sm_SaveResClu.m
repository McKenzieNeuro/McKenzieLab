function spikes =  sm_SaveResClu(basename,savepath,spk,clu,M)



xml  = LoadXml([savepath filesep basename '.xml']);

res = spk*xml.SampleRate;
% get shank assignments

[a,b] = min(M,[],2);
[~,b1] = min(a,[],1);
b1 = squeeze(b1);


ch = {xml.AnatGrps.Channels};
ch1 = cell2mat(ch);
goodCh = find(~cell2mat({xml.AnatGrps.Skip}));
goodCh  = ch1(goodCh);

shank = nan(length(b1),1);
for i = 1:length(b1)
    
    
    shank(i) = find(cellfun(@(a) any(a==b1(i)),ch));
    
end

uclu = unique(clu);
for tgroup = 1:max(shank)
    
    inclu = uclu(shank==tgroup);
    
    if any(inclu)
        kp = ismember(clu,inclu);
        tclu = clu(kp);
        tspktimes = res(kp);
        
        cluname = fullfile(savepath, [basename '.clu.' num2str(tgroup)]);
        resname = fullfile(savepath, [basename '.res.' num2str(tgroup)]);
        
        
        %clu
        % if ~exist(cluname,'file')
        tclu = cat(1,length(unique(tclu)),double(tclu(:)));
        fid=fopen(cluname,'w');
        %     fprintf(fid,'%d\n',clu);
        fprintf(fid,'%.0f\n',tclu);
        fclose(fid);
        clear fid
        % end
        
        %res
        fid=fopen(resname,'w');
        fprintf(fid,'%.0f\n',tspktimes(:));
        fclose(fid);
        clear fid
        
        
    end
end

uclu = unique(clu);
for i = 0:length(uclu)-1
spikes.times{i+1} = spk(clu==i);
end


[shank,b] = sort(shank);
spikes.shank = shank
spikes.times = spikes.times(b);



end