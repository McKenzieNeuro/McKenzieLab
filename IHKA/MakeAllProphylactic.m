% get all prophalactic stim

% make sure that we have an annotation file for each dat file

fils = getAllExtFiles('R:\DGregg\NeuralData\EDS','dat',1);
kp = contains(fils,'Prophylactic') & contains(fils,'amplifier.dat');
fils = fils(kp);

dirN = fileparts(fils);
%%
for i = 1:length(dirN)
    
    
    % get stim level, whether LDT was stimulated prior, subject name, sz
    % level
    
    sav = false;
    cd(dirN{i})
    load([fileparts(dirN{i}) filesep 'recInfo.mat'])
    for j = 1:length(recInfo.fileData)
        if ~isempty(recInfo.fileData(j).blockDuration)
            if length(recInfo.fileData(j).blockDuration{1,2})==1
                sav = true;
                for k = 1:size(recInfo.fileData(j).blockDuration,1)
                    
                    if k ==1
                        recInfo.fileData(j).blockDuration{k,2} = [1 recInfo.fileData(j).blockDuration{k,2}];
                    else
                        recInfo.fileData(j).blockDuration{k,2} = [1+recInfo.fileData(j).blockDuration{k-1,2}(2) recInfo.fileData(j).blockDuration{k-1,2}(2)+ recInfo.fileData(j).blockDuration{k,2}];
                    end
                end
                
            end
        end
    end
    
    if sav
        save([fileparts(dirN{i}) filesep 'recInfo.mat'],'recInfo','-v7.3');
    end
end
%%

allsubj = [];

uSub = [ ...
    {'EDS 1.1' } ; ...
    {'EDS 1.3' }; ...
    {'EDS 2.1' }; ...
    {'EDS 2.3' }; ...
    {'LDT 10.0'}; ...
    ];

clear rat
rat(1).name = 'EDS 1.1';
rat(2).name = 'EDS 1.3';
rat(3).name = 'EDS 2.1';
rat(4).name = 'EDS 2.3';
rat(5).name = 'LDT 10.0';
totRec = zeros(5,1);
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
            
            if size(recInfo.fileData(j).blockDuration,1)>4
                
                % define each block
                
                rat(b).ses( totRec(b)+1).block(1).sz_tim = [0 0];
                rat(b).ses( totRec(b)+1).block(2).sz_tim = [0 0];
                % get control block
                
                ep_c = cell2mat(recInfo.fileData(j).blockDuration(contains(recInfo.fileData(j).blockDuration(:,1),'SX'),2));
                ep_s = cell2mat(recInfo.fileData(j).blockDuration(contains(recInfo.fileData(j).blockDuration(:,1),'S1'),2));
                
                eps = sortrows([ep_c;ep_s]);
                ix_s = 1;
                if ~isempty(ep_s) && ~isempty(ep_c) && ep_s(1)> ep_c(1)
                    ix_s = 2;
                elseif isempty(ep_s)
                    ix_s =3;
                elseif isempty(ep_c)
                    ix_s =4;
                end
                
                %get each block
                
                if ix_s==2
                    rat(b).ses( totRec(b)+1).block(1).LDTuA = 0;
                    rat(b).ses( totRec(b)+1).block(2).LDTuA  = recInfo.fileData(j).stim1_uA;
                elseif ix_s==1
                    rat(b).ses( totRec(b)+1).block(1).LDTuA  = recInfo.fileData(j).stim1_uA;
                    rat(b).ses( totRec(b)+1).block(2).LDTuA = 0;
                elseif ix_s==3
                    rat(b).ses( totRec(b)+1).block(1).LDTuA  =0;
                    rat(b).ses( totRec(b)+1).block(2).LDTuA = 0;
                else
                    rat(b).ses( totRec(b)+1).block(1).LDTuA  =  recInfo.fileData(j).stim1_uA;
                    rat(b).ses( totRec(b)+1).block(2).LDTuA =  recInfo.fileData(j).stim1_uA;
                end
                
                rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo.fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo.fileData(j).stim2_chan(1));
                
                rat(b).ses( totRec(b)+1).block(2).inductionAmp = recInfo.fileData(j).stim2_uA;
                rat(b).ses( totRec(b)+1).block(2).inductionReg = chName(recInfo.fileData(j).stim2_chan(1));
            else
                rat(b).ses( totRec(b)+1).block(1).sz_tim = [0 0];
                ep_c = cell2mat(recInfo.fileData(j).blockDuration(contains(recInfo.fileData(j).blockDuration(:,1),'SX'),2));
                ep_s = cell2mat(recInfo.fileData(j).blockDuration(contains(recInfo.fileData(j).blockDuration(:,1),'S1'),2));
                
                eps = sortrows([ep_c;ep_s]);
                if length(recInfo.fileData(j).stim1_timeStamps) ==0
                    rat(b).ses( totRec(b)+1).block(1).LDTuA = 0;
                    rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo.fileData(j).stim2_uA;
                    rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo.fileData(j).stim2_chan(1));
                    
                else
                    rat(b).ses( totRec(b)+1).block(1).LDTuA = recInfo.fileData(j).stim1_uA;
                    rat(b).ses( totRec(b)+1).block(1).inductionAmp = recInfo.fileData(j).stim2_uA;
                    rat(b).ses( totRec(b)+1).block(1).inductionReg = chName(recInfo.fileData(j).stim2_chan(1));
                end
                
                
            end
            for k = 1:size(eps,1)
                rat(b).ses( totRec(b)+1).block(k).time = eps(k,:);
            end
            
            
            % get seizure info
            ev = LoadEvents('amplifier.evt.szr');
            
            des = cellfun(@lower,ev.description,'UniformOutput',false);
            des = cellfun(@(a) strrep(a,' ',''),des,'uni',0);
            kp_on = contains(des,su) & contains(des,'on');
            kp_off = contains(des,su) & contains(des,'off');
            sz = [ev.time(kp_on) ev.time(kp_off)];
            
            [~,block_ix] = InIntervals(sz(:,1),[eps(:,2) eps(:,2)+60]);
            for k = 1:length(block_ix)
                rat(b).ses( totRec(b)+1).block(block_ix(k)).sz_tim = sz(k,:);
            end
            rat(b).ses( totRec(b)+1).ep = recInfo.fileData(j).blockDuration;
            totRec(b) =  totRec(b)+1;
            
        end
    end
end
%%


for i = 1:length(rat)
    
    for j = 1:length(rat(i).ses)
      %  if  rat(i).ses(j).block(1).sz_tim(1)>0
           
            %load recInfo
            load([fileparts(rat(i).ses(j).session) filesep 'recInfo.mat'])
            
            
            for k = 1:2
                a = ismember(recInfo.fileData(k).subject,rat(i).name);
                
                if a
                    ix = k;
                    break
                end
                    
                    
            end
            %load seizure
            
            dur = diff(rat(i).ses(j).block(1).sz_tim,[],2);
            
            if k ==1
                
                ch = 1:length(recInfo.fileData(1).data);
                  CA1 = find(contains({recInfo.fileData(1).data.site},'CA1(R)'));
                  
            else
                ch = length(recInfo.fileData(1).data) + (1:length(recInfo.fileData(2).data));
                  CA1 = length(recInfo.fileData(1).data) +find(contains({recInfo.fileData(2).data.site},'CA1(R)'));
            end
            
            nCh = length(recInfo.fileData(1).data)+length(recInfo.fileData(2).data);
            
         
            
          
               dat  = LoadBinary([rat(i).ses(j).session filesep 'amplifier.dat'],'frequency',20000,'nchannels',nCh,...
                'channels',CA1,'start',4200-3600,'duration',3600);
            y = resample(double(dat),1250,20000);
            
            %get stims
            
            stim_ts = find(diff([0;y<-1e4])>0);
            stim_ts = stim_ts(diff([0;stim_ts])>1250);
            
            if any(stim_ts)
                idx=  repmat(-125:1250,length(stim_ts),1) + repmat(stim_ts,1,length(-125:1250));
            else
                idx =[];
            end
            % k  = gaussian2Dfilter([10000 1],1250);
            %             theta = InstAmplitude(BandpassFilter(y,1250,[5 12]));
            %             delta = InstAmplitude(BandpassFilter(y,1250,[1 4]));
            %             theta =nanconvn(theta,k);
            %             delta =nanconvn(delta,k);
            %             TD = (theta./delta);
            ok = abs(awt_freqlist(y,1250,logspace(log10(1),log10(300),100)))';
            ok(:,idx) = nan;
            ts = (1:size(ok,2))/1250;
            bins = sort(3600 - [(1./(2.^(0:11))*3600) 0]);
            
            for tt = 1:100
            sp(tt,:) = avghist(ts,ok(tt,:),bins);
            end
             rat(i).ses(j).block(1).spectra = sp(:,1:end-1);
             j
%            h =  figure;
%             for dd = 1:size(dat,2)
%                 
%                 plot(double(dat(:,dd))-dd*6000,'k')
%                 hold on
%             end
%             
%              x = input('gen');
%              x  =str2num(x);
%              if x==1
%                  rat(i).ses(j).block(1).gen = 1;
%              else
%                   rat(i).ses(j).block(1).gen = 0;
%              end
%            title(rat(i).name)
%            print(h, '-dpsc2','E:\Dropbox\UNM\Presentations\2023\ParkCity\evokesz.ps' ,'-append')
%          
%             close all
      %  end
    end
i    
end
%sm_ps2pdf('E:\Dropbox\UNM\Presentations\2023\ParkCity\evokesz.ps','E:\Dropbox\UNM\Presentations\2023\ParkCity\evokesz.pdf')

%%
load('R:\DGregg\NeuralData\EDS\prophylacticStim.mat')

clear tstat
for f = 1:100
for t = 1:12
    ix_s = ones(5,1);
    ix_ns = ones(5,1);
    clear nostim stim inductAmp_s inductAmp_ns gen_s gen_ns  TD_s TD_ns
    for i = 1:5
        
        for k = 1:length(rat(i).ses)
            if ~isfield(rat(i).ses(k).block(1),'gen') 
                rat(i).ses(k).block(1).gen = 0;
            end
            if rat(i).ses(k).block(1).LDTuA>0 && ~strcmp(rat(i).ses(k).block(1).inductionReg{1},'AVT(L)')% &&  rat(i).ses(k).block(1).inductionAmp==400
                stim{i,ix_s(i)} =  rat(i).ses(k).block(1).sz_tim;
                inductAmp_s(i,ix_s(i)) = rat(i).ses(k).block(1).inductionAmp;
                gen_s{i,ix_s(i)} =  rat(i).ses(k).block(1).gen;
               % TD_s{i,ix_s(i)}  = rat(i).ses(k).block(1).TD(t);
                TD_s{i,ix_s(i)}  =   rat(i).ses(k).block(1).spectra(f,t);
              
                
                ix_s(i) = ix_s(i)+1;
            elseif  rat(i).ses(k).block(1).LDTuA==0 && ~strcmp(rat(i).ses(k).block(1).inductionReg{1},'AVT(L)')% &&  rat(i).ses(k).block(1).inductionAmp==400
                
                nostim{i,ix_ns(i)} =  rat(i).ses(k).block(1).sz_tim;
                
                gen_ns{i,ix_ns(i)} =  rat(i).ses(k).block(1).gen;
             %   TD_ns{i,ix_ns(i)}  = rat(i).ses(k).block(1).TD(t);
                 TD_ns{i,ix_ns(i)}  =   rat(i).ses(k).block(1).spectra(f,t);
                ix_ns(i) = ix_ns(i)+1;
            end
        end
    end
    
    
    
    
    stim = cell2mat(cellfun(@(a) nanmean(diff(a,[],2)),stim,'UniformOutput',false));
    nostim = cell2mat(cellfun(@(a) nanmedian(diff(a,[],2)),nostim,'UniformOutput',false));
    
    
    stim_gen = cellfun(@(a) nanmean(a),gen_s);
    nostim_gen = cellfun(@(a) nanmean(a),gen_ns);
    
    stim_TD = cellfun(@(a) nanmean(a),TD_s);
    nostim_TD = cellfun(@(a) nanmean(a),TD_ns);
    
    
    
    subjID1 = repmat([1:5]',1,size(stim,2));
    subjID2 = repmat([1:5]',1,size(nostim,2));
    gr1 = ones(size(stim));
    gr2 = 2*ones(size(nostim));
    
    norm_stim = stim./max([stim nostim],[],2);
    norm_nostim = nostim./max([stim nostim],[],2);
    gen = [stim_gen(:);nostim_gen(:)];
    sz_dur = [stim(:);nostim(:)];
    TD = [stim_TD(:);nostim_TD(:)];
    gr = [gr1(:);gr2(:)];
    subjID = [subjID1(:);subjID2(:)];
   % TD_all(:,t) = TD;
    tbl = table(sz_dur,TD,gen,gr,subjID,'VariableNames',{'seiz_dur','TD','gen','stim','subjID'});
    kp= ~isnan(TD) & ~isnan(gen);
    lme = fitglme(tbl(kp,:),'gen ~ TD + (1+TD|subjID)','distribution','binomial');
    tstat(f,t,:) = [lme.Coefficients.tStat(2) lme.Coefficients.pValue(2)];
    t
end
end
%%

for i = 1:5
    maxTD(i) = max(max(TD_all(subjID==i,:)));
    
    TD_all_norm(subjID==i,:)=  TD_all(subjID==i,:)/maxTD;
end

%%
clear tstat
for i = 1:size(TD_all_norm,2)
 tbl = table(sz_dur,TD_all_norm(:,i),gen,gr,subjID,'VariableNames',{'seiz_dur','TD','gen','stim','subjID'});
   lme = fitglme(tbl(kp,:),'gen ~ TD + (1+TD|subjID)','distribution','binomial');
    tstat(i,:) = [lme.Coefficients.tStat(2) lme.Coefficients.pValue(2)];
end

figure

semilogx(3600-bins(1:end-1),nanmean(TD_all(gen==1,:)),'w')
hold on
plotMeanSEM(3600-bins(1:end-1),TD_all_norm(gen==0,:),'r')
plotMeanSEM(3600-bins(1:end-1),TD_all_norm(gen==1,:),'k')
ylim([0 .6])
   
%%
% check each seizure 
semilogx(3600-bins(1:end-1),(tstat(:,1)))


%%
close all
col = linspecer(5,'jet');
for i = 1:5
    
    plot(gr(subjID==i)+(rand(sum(subjID==i),1)-.5)/10,sz_dur(subjID==i),'.','markersize',20,'color',col{i})
    hold on
end
figure
for i = 1:5
    
    plot(gr(subjID==i)+(rand(sum(subjID==i),1)-.5)/10,TD(subjID==i),'.','markersize',20,'color',col{i})
    hold on
end




figure
for i = 1:5
plot([1 2],[nanmean(stim_gen(i,:),2) nanmean(nostim_gen(i,:),2)],'color',col{i})
hold on
end

%%

%plot example
cd('R:\DGregg\NeuralData\EDS\Prophylactic3\4-14-2023(12.59)\RHS_230414_130000')

d = LoadBinary('amplifier.dat','nchannels',15,'channels',1:7,'frequency',20000,'start',8397,'duration',100);

%%

figure
for i = 1:7
    y = resample(double(d(:,i)),1250,20000);
    ts = (1:length(y))/1250-(8411.4-8397);
    plot(ts,y*.195 - i*1000,'k')
    hold on
end

%%
close all
figure

ts = sort(-(1:length((1250*6:1250*14.5)))/1250)
plot(ts,.195*y(1250*6:1250*14.5))
hold on
plot(ts,.195*BandpassFilter(y(1250*6:1250*14.5),1250,[5 12]),'r')
ylim([-400 400])

    %%
    
cd('R:\DGregg\NeuralData\EDS\Prophylactic6\5-2-2023(12.58)\RHS_230502_130000')
    
d = LoadBinary('amplifier.dat','nchannels',16,'channels',1:8,'frequency',20000,'start',4201-14.4,'duration',100);



figure
for i = 1:8
    y = resample(double(d(:,i)),1250,20000);
    ts = (1:length(y))/1250-14.4;
    plot(ts,y*.195 - i*1000,'k')
    hold on
end

%%
close all
figure

ts = sort(-(1:length((1250*5.5:1250*14)))/1250)
plot(ts,.195*y(1250*5.5:1250*14))
hold on
plot(ts,.195*BandpassFilter(y(1250*5.5:1250*14),1250,[5 12]),'r')
ylim([-400 400])