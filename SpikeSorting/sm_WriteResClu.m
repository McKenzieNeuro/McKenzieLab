function sm_WriteResClu(spikes,basepath)

basename = spikes.sessionName;
groups = spikes.shankID;

for tgroup = unique(groups)
%%
cluname = fullfile(basepath, [basename '.clu.' num2str(tgroup)]);
ix = spikes.shankID==tgroup;
tclu = spikes.UID(ix);
tclu = cellfun(@(a,b) a*ones(length(b),1),num2cell(tclu),spikes.times(ix),'uni',0);
tclu = cell2mat(tclu');
tspktimes = round(cell2mat(spikes.times(ix)')*30000);
[tspktimes,b] = sort(tspktimes);
tclu = tclu(b);

tclu = cat(1,length(unique(tclu)),double(tclu));
fid=fopen(cluname,'w');
%     fprintf(fid,'%d\n',clu);
fprintf(fid,'%.0f\n',tclu);
fclose(fid);
clear fid
%%

%res


resname = fullfile(basepath, [basename '.res.' num2str(tgroup)]);
fid=fopen(resname,'w');
fprintf(fid,'%.0f\n',tspktimes);
fclose(fid);
clear fid
end