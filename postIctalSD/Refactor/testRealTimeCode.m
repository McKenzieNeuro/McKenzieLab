    model = load('R:\Analysis\McKenzieLab\postIctalSD\Refactor\Models\model_pilo_bilateral_24FEB26_sk.mat');
    dSF=1250;    
    nChannel = length(model.ops.ch_subj);
    freq = model.ops.freqs;
    nFreq = length(freq);
    feat = nan(1,nFreq*nChannel);

    modelOPs.freqs = freq;
    modelOPs.durFeat = model.ops.durFeat; % 4s feature bins
    modelOPs.art_thres = 5e4; % NOTE probably needs to be decreased
    modelOPs.features = model.ops.features;
    modelOPs.ch_subj = model.ops.ch_subj;
    modelOPs.Fs = dSF;


while 1                    
    if spmdProbe(intanIO,ampData) == 1

        %data = spmdReceive(intanIO,ampData); 
        path="E:\Recordings\v8_6\2-25-2026(11.7)\2-25-2026(11.14)\RHS_260225_111605\amplifier.dat"
        data = LoadBinary(path,'nchannels',7,'frequency', ...
                20000,'channels',[1,5],'duration', 4,'start',0);
        if ~isempty(data)

            data = data'

            disp(size(data))
            sampleRatio=sF/modelOPs.Fs;
            ratio = 450/(sF/2) ;

            %ntbuff should be even multiple of sampleRatio
            ntbuff = 525;  %default filter size in iosr toolbox
            disp('ntbuff');
            if mod(ntbuff,sampleRatio)~=0
                ntbuff = ntbuff + sampleRatio-mod(ntbuff,sampleRatio);
            end
            
            ds_data=[[],[]];
            disp('filt tmp');
            tmp = iosr.dsp.sincFilter(data',ratio);
            disp('transpose tmp');
            tmp=tmp';
            
            for i = 1:size(data,1)
                disp(i)
                if i==1
                    disp('ds tmp pass 1')
                elseif i==2
                    disp('ds tmp pass 2')
                end
                ds_data(i,:) = tmp(i,ntbuff+sampleRatio:sampleRatio:end);
            
            end

            disp(['ntbuff: ', num2str(ntbuff), ' sampleRatio: ',...
                num2str(sampleRatio),' ratio: ',num2str(ratio)])
            %this will need to be changed if feature
            %changes
            %if ~any(any(abs(data) > model.ops.art_thres))
            if ~any(any(abs(ds_data) > inf))
                 feat = sm_GetDataFeature_rat(ds_data',[],modelOPs);

                 label = predict(model.rusTree,feat);