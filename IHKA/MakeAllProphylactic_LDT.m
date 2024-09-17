%ipsi (red0)
dirN = [ ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(12.58)\RHS_P2_231202_150025'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-4-2023(12.58)\RHS_P2_231204_150025'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(12.58)\RHS_P2_231202_150025'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(12.58)\RHS_P2_231202_150025'} ; ...
    

% %contra LFS

    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-21-2023(16.53)\RHS_231121_165356'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-22-2023(12.58)\RHS_231122_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-23-2023(12.58)\RHS_231123_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(12.58)\RHS_231124_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-25-2023(12.58)\RHS_231125_130000'} ; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-27-2023(13.4)\RHS_231127_130449'} ; ...
    
    
    

{'R:\DGregg\NeuralData\EDS_Cohort2\Induction3\11-12-2023(12.58)\RHS_231112_130000'  } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction3\11-10-2023(13.7)\RHS_231110_130745'  } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction2\11-9-2023(12.58)\RHS_231109_130000'  } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction2\11-8-2023(13.41)\RHS_231108_134150'  } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction2\11-7-2023(13.31)\RHS_231107_133223'  } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction1\11-4-2023(12.58)\RHS_231104_130000'  } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction1\11-3-2023(14.26)\RHS_231103_142651' } ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction1\11-2-2023(12.58)\RHS_231102_130000'} ; ...
]
%%


allsubj = [];

uSub = [ ...
    {'EDS 4.0' } ; ...
    {'EDS 2.3' }; ...
    {'EDS 5.1' }; ...
    {'EDS 3.0' }; ...
     {'EDS 4.1' }; ...
       {'EDS 4.2' }; ...
    ];

clear rat
rat(1).name = 'EDS 4.0';
rat(2).name = 'EDS 2.3';
rat(3).name = 'EDS 5.1';
rat(4).name = 'EDS 3.0';
rat(5).name = 'EDS 4.1';
rat(6).name = 'EDS 4.2';

totRec = zeros(6,1);
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
        
        if iscell(recInfo(1).fileData(j).subject)
            tmp = recInfo(1).fileData(j).subject{1};
            
            
            subjs = [subjs;num2cell(tmp,2)];
        end
    end
    
    for j = 1:length(subjs)
        
        if ~all(contains(recInfo(1).fileData(j).blockDuration(:,1),'BL')) && ~all(contains(recInfo(1).fileData(j).blockDuration(:,1),'XX'))
            
            [~,b] = ismember(recInfo(1).fileData(j).subject,uSub);
            
            
            rat(b).ses( totRec(b)+1).session = dirN{i};
            su = recInfo(1).fileData(j).subject{1};
            su = lower(strrep(su,' ',''));
            %one or two block
            
            
            rat(b).ses( totRec(b)+1).block(1).sz_tim = [0 0];
            ep_c = cell2mat(recInfo(1).fileData(j).blockDuration(contains(recInfo(1).fileData(j).blockDuration(:,1),'S1_O'),2));
            ep_s = cell2mat(recInfo(1).fileData(j).blockDuration(contains(recInfo(1).fileData(j).blockDuration(:,1),'S2'),2));
            
            eps = sortrows([ep_c;ep_s]);
            if ~isfield(recInfo(1).fileData(j),'stim1_timeStamps')
                
                if contains(recInfo(1).fileData(j).blockDuration(:,1),'S1_O')
                    rat(b).ses( totRec(b)+1).block(1).LDTuA = recInfo(1).fileData(j).stim1_uA;
                else
                    rat(b).ses( totRec(b)+1).block(1).LDTuA =0;
                end
                rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo(1).fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo(1).fileData(j).stim2_chan(1));
                
             elseif length(recInfo(1).fileData(j).stim1_timeStamps) ==0
                rat(b).ses( totRec(b)+1).block(1).LDTuA = 0;
                rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo(1).fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo(1).fileData(j).stim2_chan(1));
                
            else
                rat(b).ses( totRec(b)+1).block(1).LDTuA = recInfo(1).fileData(j).stim1_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo(1).fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo(1).fileData(j).stim2_chan(1));
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
            
            rat(b).ses( totRec(b)+1).ep = recInfo(1).fileData(j).blockDuration;
            
            tp = any(contains(recInfo(1).fileData(1).blockDuration(:,1),'S2_ptC'));
            rat(b).ses( totRec(b)+1).block.thetaPeak = tp;
            totRec(b) =  totRec(b)+1;
            
        end
    end
end
%%
sz_len = cell(6,4);
sz_gen = cell(6,4);
for i = 1:length(rat)
    
    for j = 1:length(rat(i).ses)
        
        if rat(i).ses(j).block.LDTuA==0 & rat(i).ses(j).block.thetaPeak==0
            cond = 1;
        elseif  rat(i).ses(j).block.LDTuA>0 & rat(i).ses(j).block.thetaPeak==0
            cond = 2;
        elseif  rat(i).ses(j).block.LDTuA==0 & rat(i).ses(j).block.thetaPeak>0
            cond = 3;
        elseif rat(i).ses(j).block.LDTuA>0 & rat(i).ses(j).block.thetaPeak>0
            cond = 4;
        end
        
        tmp  = diff(rat(i).ses(j).block.sz_tim,[],2);
        
        if isfield(rat(i).ses(j).block,'gen')rat(i).ses(j).block
        gen = rat(i).ses(j).block.gen;
        else
            gen = nan;
        end
        if tmp<5
            tmp = nan;
            gen = nan;
        end
         sz_len{i,cond} = [sz_len{i,cond};tmp];
         
          sz_gen{i,cond} = [sz_gen{i,cond};gen];
    end
  
   
end
%%
close all

for i = 6:length(rat)
    
    for j = 1:length(rat(i).ses)
      %  if  rat(i).ses(j).block(1).sz_tim(1)>0
           
            %load recInfo
            load([fileparts(rat(i).ses(j).session) filesep 'recInfo.mat'])
            
              
           
            %load seizure
            
            dur = diff(rat(i).ses(j).block(1).sz_tim,[],2);
            for k = 1:4
                
                if iscell(recInfo(1).fileData(k).subject)
                     a = ismember(recInfo(1).fileData(k).subject{1},rat(i).name);
                else
                    
                    a = ismember(recInfo(1).fileData(k).subject,rat(i).name);
                end
                if a
                    ix = k;
                    break
                end
                    
                    
              end
              
               ch_all =  {recInfo(1).fileData.data};
               all_prev = cumsum(cellfun(@length,ch_all));
                 ch = 1:length(ch_all{k});
           if k >1
               ch = ch+all_prev(k-1);
           end
              
            nCh = sum(cellfun(@length,{recInfo(1).fileData.data}));
            
         
            
          if ~isempty(rat(i).ses(j).block(1).sz_tim)
               dat  = LoadBinary([rat(i).ses(j).session filesep 'amplifier.dat'],'frequency',20000,'nchannels',nCh,...
                'channels',ch,'start',rat(i).ses(j).block(1).sz_tim(1)-5,'duration',100);
           % y = resample(double(dat),1250,20000);
            
            %get stims
            
            %stim_ts = find(diff([0;y<-1e4])>0);
            %stim_ts = stim_ts(diff([0;stim_ts])>1250);
            
          %  if any(stim_ts)
          %      idx=  repmat(-125:1250,length(stim_ts),1) + repmat(stim_ts,1,length(-125:1250));
           % else
           %     idx =[];
           % end
            % k  = gaussian2Dfilter([10000 1],1250);
            %             theta = InstAmplitude(BandpassFilter(y,1250,[5 12]));
            %             delta = InstAmplitude(BandpassFilter(y,1250,[1 4]));
            %             theta =nanconvn(theta,k);
            %             delta =nanconvn(delta,k);
            %             TD = (theta./delta);
          %  ok = abs(awt_freqlist(y,1250,logspace(log10(1),log10(300),100)))';
          %  ok(:,idx) = nan;
          %  ts = (1:size(ok,2))/1250;
          % bins = sort(3600 - [(1./(2.^(0:11))*3600) 0]);
            
           % for tt = 1:100
           % sp(tt,:) = avghist(ts,ok(tt,:),bins);
           %end
            % rat(i).ses(j).block(1).spectra = sp(:,1:end-1);
             
            
            
           
%            h =  figure;
            for dd = 1:size(dat,2)
                
                plot(double(dat(:,dd))-dd*6000,'k')
                hold on
            end
            
             x = input('gen');
           %  x  =str2num(x);
             if x==1
                 rat(i).ses(j).block(1).gen = 1;
             else
                  rat(i).ses(j).block(1).gen = 0;
             end
%            title(rat(i).name)
%            print(h, '-dpsc2','E:\Dropbox\UNM\Presentations\2023\ParkCity\evokesz.ps' ,'-append')
%          
            close all
      %  end
          end
    end
i    
end

save('R:\DGregg\NeuralData\EDS_Cohort2\all_theta_stim.mat','rat')
%%
kp = all(cellfun(@(a) all(a==0),sz_len),2);
sz_len(kp,:) = [];

