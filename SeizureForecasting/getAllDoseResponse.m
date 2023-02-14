% thalamus
topDirThal = 'R:\DGregg\NeuralData\LTD 7.0\9-1-2022(15.20)_AVT(L)_RHS';
order1= [{'PrL(L)'}    {'PrL(R)'}    {'AVT(L)'}    {'CA1(L)'}    {'CA1(R)'}];
%CA1
topDirCA1 = 'R:\DGregg\NeuralData\LTD 7.0\8-23-2022(14.7)_CA1(L)_RHS';
order1 = [{'PrL(L)'}    {'PrL(R)'}    {'AVT(L)'}    {'CA1(L)'}    {'CA1(R)'}];
% BLA
topDirBLA = 'R:\DGregg\NeuralData\LTD 9.0\8-22-2022(13.38)_BLA(R)_RHS';
order2 = [ {'PrL(L)'}    {'PrL(R)'}    {'AVT(L)'}    {'BLA(R)'}    {'CA1(R)'}    {'gRSC(L)'}    {'LTD(L)'}];
% LTD
topDirLTD = 'R:\DGregg\NeuralData\LTD 9.0\8-18-2022(16%48)_LTD(L)_RHS';
order2 = [ {'PrL(L)'}    {'PrL(R)'}    {'AVT(L)'}    {'BLA(R)'}    {'CA1(R)'}    {'gRSC(L)'}    {'LTD(L)'}];

allReg = [{'PrL(L)'}     {'AVT(L)'}   {'gRSC(L)'}  {'CA1(L)'}    {'LTD(L)'} {'PrL(R)'}  {'CA1(R)'}  {'BLA(R)'}  ];

%%
%clear d
for i = 2
    
    switch i
        case 1
            topDir = topDirThal;
        case 2
            topDir = topDirCA1;
        case 3
            topDir = topDirLTD;
        case 4
            topDir = topDirBLA;
            
    end
    cd(topDir)
    ok = dir(topDir);
    subDirs = {ok(contains({ok.name},'uA') &  cell2mat({ok.isdir})).name}';
    
    for j = 1:length(subDirs)
        
        cd(subDirs{j})
        if exist('TTL_pulse.mat')
            v=load('TTL_pulse.mat');
            TTL = v.ups;
        else
            TTL = sm_getDigitalin(pwd,'digitalin',1,20000);
        end
        out = read_Intan_RHS2000_file('info.rhs');
        nCh = out.num_amplifier_channels;
        
        dt = nan(nCh,40000,length(TTL));
        for k = 1:length(TTL)
            dt(:,:,k) = LoadBinary('amplifier.dat','nchannels',nCh,'channels',1:nCh,'frequency',20000,'start',TTL(k)-1,'duration',2)';
            
        end
        
        d{i,j} = nanmean(dt,3);
        cd(topDir)
        [i j]
    end
end





%%
close all
pow  =3;
figure
ax = tight_subplot(4,8);

for i = 1:numel(ax)
    
    axes(ax(i))
    axis off
end

for i = 1:4
    
    if i<3
        for j = 1:5
            tmp = d{i,pow}(j,:);
            tmp(20000:20750) = nan;
            [a,b]  =ismember(order1,allReg);
            axes(ax(i,b(j)))
            axis on
            plot((1:length(d{i,pow}(j,:)))/20000 -1 , 0.195 *tmp)
            ylim([-5000 5000])
            xlim([-.25 1])
            
        end
    else
        for j = 1:7
            tmp = d{i,pow}(j,:);
            tmp(20000:20750) = nan;
            [a,b]  =ismember(order2,allReg);
            axes(ax(i,b(j)))
              axis on
            plot((1:length(d{i,pow}(j,:)))/20000 -1 , 0.195 *tmp)
            xlim([-.25 1])
            ylim([-4000 4000])
            
        end
    end
    
    
    
    
    
end

for i = 1:8
    axes(ax(1,i))
    title(allReg(i))
    for j = 1:4
         axes(ax(j,i))
    if i >1
        set(gca,'yticklabel',[])
    end
    end
end



for i = 1:4
    axes(ax(i,1))
    switch i
        case 1
            lab =  'stim AVT (L)';
        case 2
            lab =  'stim CA1 (L)';
        case 3
            lab =  'stim LDT (L)';
        case 4
            lab =  'stim BLA (R)';
    end
    ylabel(lab)
    
   for j = 1:8
       axes(ax(i,j))
    if i < 4
             set(gca,'xticklabel',[])
    end
   end
end

%%
close all
ax = tight_subplot(1,5)
for j = 1:5
    axes(ax(j))

  col = linspecer(9,'jet');
    
for i = 2:9
  
    plot((1:length(d{1,i}(1,:)))/20000 -1,.195*d{2,i}(j,:),'color',col{i},'linewidth',2)
    hold on
     xlim([-.25 1])
            ylim([-4000 4000])
end


colormap(jet(9))
cbh = colorbar ; %Create Colorbar
 cbh.Ticks = linspace(0, 1, 9) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = num2cell(0:25:200) ;  

xlabel('time from CA1 stim (s)')
ylabel('response (mV)')
ylabel(cbh,'Stim Level(uA)')
title(order1{j})
end
