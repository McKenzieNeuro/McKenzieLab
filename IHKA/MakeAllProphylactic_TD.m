%ipsi (red0)
% LFSfils = [ ...
%     {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(12.58)\RHS_P2_231202_150025'} ; ...
%     {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-4-2023(12.58)\RHS_P2_231204_150025'} ; ...
%     {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(12.58)\RHS_P2_231202_150025'} ; ...
%     {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(12.58)\RHS_P2_231202_150025'} ; ...
%     ];

%%

dirN = [ ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-21-2023(16.53)\RHS_231121_165356'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-22-2023(12.58)\RHS_231122_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-23-2023(12.58)\RHS_231123_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(12.58)\RHS_231124_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-25-2023(12.58)\RHS_231125_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-27-2023(13.4)\RHS_231127_130449'} ; ...
    
    
    
    ];

%%


allsubj = [];

uSub = [ ...
    {'EDS 4.0' } ; ...
    {'EDS 2.3' }; ...
    {'EDS 5.1' }; ...
    {'EDS 3.0' }; ...
    
    ];

clear rat
rat(1).name = 'EDS 4.0';
rat(2).name = 'EDS 2.3';
rat(3).name = 'EDS 5.1';
rat(4).name = 'EDS 3.0';


totRec = zeros(4,1);
chName = [...
    {'M2(L)' } ; ...
    {'M2(R)' } ; ...
    {'AVT(L)' } ; ...
    {'BLA(R)' } ; ...
    {'CA1(L)' } ; ...
    {'CA1(R)' } ; ...
    {'LDT(L1)'} ; ...
    {'LDT(L2)'} ; ...
    ];
warning off
for i = 1:length(dirN)
    
    
    % get stim level, whether LDT was stimulated prior, subject name, sz
    % level
    
    
    cd(dirN{i})
    load([fileparts(dirN{i}) filesep 'recInfo.mat'])
    %get subject names
    subjs =[];
    for j = 1:4
        
        if iscell(recInfo.fileData(j).subject)
            tmp = recInfo.fileData(j).subject{1};
            
            
            subjs = [subjs;num2cell(tmp,2)];
        end
    end
    
    for j = 1:length(subjs)
        
        if ~all(contains(recInfo.fileData(j).blockDuration(:,1),'BL')) && ~all(contains(recInfo.fileData(j).blockDuration(:,1),'XX'))
            
            [~,b] = ismember(recInfo.fileData(j).subject,uSub);
            
            
            rat(b).ses( totRec(b)+1).session = dirN{i};
            su = recInfo.fileData(j).subject{1};
            su = lower(strrep(su,' ',''));
            %one or two block
            
            
            rat(b).ses( totRec(b)+1).block(1).sz_tim = [0 0];
            ep_c = cell2mat(recInfo.fileData(j).blockDuration(contains(recInfo.fileData(j).blockDuration(:,1),'S1_O'),2));
            ep_s = cell2mat(recInfo.fileData(j).blockDuration(contains(recInfo.fileData(j).blockDuration(:,1),'S2'),2));
            
            eps = sortrows([ep_c;ep_s]);
            if length(recInfo.fileData(j).stim1_timeStamps) ==0
                rat(b).ses( totRec(b)+1).block(1).CA1uA = 0;
                rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo.fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo.fileData(j).stim2_chan(1));
                
            else
                rat(b).ses( totRec(b)+1).block(1).CA1uA = recInfo.fileData(j).stim1_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo.fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo.fileData(j).stim2_chan(1));
            end
            
            
            %
            %             for k = 1:size(eps,1)
            %                 rat(b).ses( totRec(b)+1).block(k).time = eps(k,:);
            %             end
            
            
            % get seizure info
            ev = LoadEvents('amplifier.evt.szr');
            
            des = cellfun(@lower,ev.description,'UniformOutput',false);
            des = cellfun(@(a) strrep(a,' ',''),des,'uni',0);
            kp_on = contains(des,su) & contains(des,'on');
            kp_off = contains(des,su) & contains(des,'off');
            sz = [ev.time(kp_on) ev.time(kp_off)];
            
            rat(b).ses( totRec(b)+1).block(1).sz_tim = sz;
            
            rat(b).ses( totRec(b)+1).ep = recInfo.fileData(j).blockDuration;
            
            tp = any(contains(recInfo.fileData(1).blockDuration(:,1),'S2_ptC'));
            rat(b).ses( totRec(b)+1).block.thetaPeak = tp;
            totRec(b) =  totRec(b)+1;
            
        end
    end
end
%%
sz_len = cell(4,4);
for i = 1:length(rat)
    
    for j = 1:length(rat(i).ses)
        
        if rat(i).ses(j).block.CA1uA==0 & rat(i).ses(j).block.thetaPeak==0
            cond = 1;
        elseif  rat(i).ses(j).block.CA1uA>0 & rat(i).ses(j).block.thetaPeak==0
            cond = 2;
        elseif  rat(i).ses(j).block.CA1uA==0 & rat(i).ses(j).block.thetaPeak>0
            cond = 3;
        elseif rat(i).ses(j).block.CA1uA>0 & rat(i).ses(j).block.thetaPeak>0
            cond = 4;
        end
        
        tmp  = diff(rat(i).ses(j).block.sz_tim,[],2);
        
        if tmp<5
            tmp = 0;
        end
         sz_len{i,cond} = [sz_len{i,cond};tmp];
    end
  
   
end

%%
kp = all(cellfun(@(a) all(a==0),sz_len),2);
sz_len(kp,:) = [];

