function [estimateLabel,trueLabel,ts1] = sm_getSeizurePred(fname,seizFil,rusTree)

TSname1 = 'Seizure starts'; % TS name in file.

TSname2 = 'ends';
TSdata = readtable(seizFil);
TSdata = table2cell(TSdata);
seizure_start = cell2mat(TSdata(cellfun(@any,regexp(TSdata(:,6),TSname1)),4));


seizure_end = cell2mat(TSdata(cellfun(@any,regexp(TSdata(:,6),TSname2)),4));



freqs = logspace(log10(.5),log10(200),20);
Fs =2000;
bins = [10 600 3600 3*3600];
powerFil = [fname '_1.dat'];
s = dir(powerFil);
dur = s.bytes/41/Fs/2;
estimateLabel =[];
dat1 = [];

ts = 1:dur;
ts1 = ts;
kp = true(size(ts));
for i = 1:length(seizure_start)
    ts1(kp&ts<seizure_start(i)) = ts(kp&ts<seizure_start(i)) - seizure_start(i);
    ts1(kp& ts>seizure_start(i) & ts<seizure_end(i)) = .5;
    ts1(kp& ts>seizure_end(i) & ts<seizure_end(i)+600) = 1.5;
    
    
    kp(ts<seizure_start(i) | ...
        (ts>seizure_start(i) & ts<seizure_end(i)) | ...
        (ts>seizure_end(i) & ts<seizure_end(i)+600)) = false;
end

[~,trueLabel] = histc(ts1,[-inf -fliplr(bins) 0 1 2]);




for i = 1:dur
    
    try
        
        datT = [];rD = nan(10000,4);
        for ch = 1:4
            powerFil = [fname '_' num2str(ch) '.dat'];
            tmp = LoadBinary(powerFil,'nchannels',41,'frequency',Fs,'channels',1:41,'duration', 5,'start',i);
            datT = [datT 1000*mean(tmp(:,2:2:39))];
            rD(:,ch) = tmp(:,1);
            
            
            if ch==2
                comod =[];
                filtered_phase = double(tmp(:,3:2:41))/1000;
                wavespec_amp = double(tmp(:,2:2:40));
                
                for apr = 15:20
                    
                    
                    
                    
                    
                    for idx = 1:10
                        
                        
                        comod(apr-14,idx) = circ_corrcl(filtered_phase(:,idx),wavespec_amp(:,apr));
                        
                    end
                    
                end
            end
        end
        cxy = nan(20,6); ix=1;
        for ch1 = 1:4
            for ch2 = ch1+1:4
                [cxy(:,ix),f] = mscohere(rD(:,ch1),rD(:,ch2),[],[],freqs,2000);
                ix= ix+1;
            end
        end
        datT = [datT comod(:)' cxy(:)'];
        dat1 = [dat1;datT];
        
        if mod(i,100)== 0
            estimateLabel = [estimateLabel;predict(rusTree,dat1)];
            dat1 =[];
        elseif i > (dur- mod(dur,100))
            estimateLabel = [estimateLabel;predict(rusTree,dat1)];
            
        end
        
    end
    
    
    
    
end



end

