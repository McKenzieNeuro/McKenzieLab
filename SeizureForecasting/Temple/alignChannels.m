topDir = 'G:\data\isip\oneTreeData\binaries';

fils = getAllExtFiles(topDir,'txt',0);

uFils = unique(cellfun(@(a) a(1:end-7),fils,'uni',0));



channel  = [ ...
    {'fp1'} ; ...
    {'fp2'} ; ...
    {'f7'} ; ...
    {'f8'} ; ...
    {'t3'} ; ...
    {'t4'} ; ...
    {'t5'} ; ...
    {'t6'} ; ...
    {'o1'} ; ...
    {'o2'} ; ...
    {'f3'} ; ...
    {'f4'} ; ...
    {'c3'} ; ...
    {'c4'} ; ...
    {'p3'} ; ...
    {'p4'} ; ...
    {'a1'} ; ...
    {'a2'} ; ...
    {'fpz'} ; ...
    {'fz'} ; ...
    {'cz'} ; ...
    {'pz'} ; ...
    {'oz'} ; ...
    {'ekg'} ; ...
    
    ];




for i = 1:length(uFils)
    ch_txt = fils(contains(fils,uFils{i}));
    
    for j = 1:length(ch_txt)
        fid = fopen(ch_txt{j});
        tline = fgetl(fid);
        
        fclose(fid);
        
        labelInfo(i).lab{j}  =tline(4:end);
    end
    
    i
    
end

%%

for i = 1:length(labelInfo)
    
    
    for j = 1:24
        tmp =  find(cellfun(@any,regexpi(labelInfo(i).lab',channel{j})),1,'first');
        if ~isempty(tmp)
            chOrder(i,j) = tmp;
        else
            chOrder(i,j) = 0;
        end
    end
    
    
end

%%

[a,b,c] = unique(chOrder,'rows');
[~,b] = ismember(chOrder,a,'rows');
