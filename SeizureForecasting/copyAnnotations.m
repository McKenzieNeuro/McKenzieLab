fils = getAllExtFiles('Y:\','LOG',1);
[a,b] = fileparts(fils);

%%

for i = 22361:length(fils)
    
    copyfile(fils{i},['Z:\UNMH_EEGCorpus\annotations\' b{i} '.txt'],'f')
    
    i
end

%%
% get all local annotations
fils = getAllExtFiles('Z:\UNMH_EEGCorpus\annotations','txt',1);

for i = 1:length(fils)
    try
        hdr = sm_ReadNKLogFile(fils{i});
        time{i} = hdr.logs.time;
        label{i} = hdr.logs.label;
    end
    i
end

%%
for i = 1:length(label)
    if ~isempty(label{i})
       gdLabel(i) =  any(cellfun(@any,regexpi(label{i},'sz'))) |  any(cellfun(@any,regexpi(label{i},'seiz')));
       
       %gdLabel(i) = any(cellfun(@any,regexpi(label{i},'Sz_4+10')));
       i
    end
end


gdLabel = find(gdLabel);
%%

for i = 1:length(gdLabel)
    l = label{gdLabel(i)};
      t = time{gdLabel(i)};
      
      e410 = contains(l,'4') & contains(l,'10');
      seiz_tim{i} =[];
        seiz_lab{i} =[];
    for j = 1:length(l)
        if ~e410(j) && (~isempty(regexpi(l{j},'sz')) | ~isempty(regexpi(l{j},'seiz'))) & ( isempty(regexpi(l{j},'no')) & isempty(regexpi(l{j},'end')) & isempty(regexpi(l{j},'stop')) & isempty(regexpi(l{j},'ove')) & isempty(regexpi(l{j},'(Per')))
          
            
            seiz_tim{i} = [seiz_tim{i};t(j)];
            
          
            seiz_lab{i} = [seiz_lab{i};l(j)];
        end
    end
    
%     if any(diff(seiz_tim{i})>999 & diff(seiz_tim{i})<1001)
%         error('here')
%     end
%    i 
i
end

s = ((cell2mat(cellfun(@diff,seiz_tim,'UniformOutput',false)')));
s = s(s>60);
%%


uL = cellfun(@unique,seiz_lab','UniformOutput',false);
LL = [];

for i = 1:length(uL)
    LL = [ LL;uL{i}];
end
%%

kp = cellfun(@any,seiz_tim);
seiz_tim = seiz_tim(kp);
