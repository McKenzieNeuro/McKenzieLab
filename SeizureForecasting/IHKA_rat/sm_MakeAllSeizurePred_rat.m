% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutputDir =  'R:\Analysis\SeizureForecasting\IHKA_rat_RF\Classification\';



FeatureFileOutput = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\features.mat';
load(FeatureFileOutput)
load('R:\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat')
ops.subjName = 'EDS4.0';
%%

warning off
% loop over files to predict
for i = 1:size(sessions,1)
   % ClassifierFileOutput = [ClassifierFileOutputDir filesep 'classification_' num2str(i) '.mat'];
   % load(ClassifierFileOutput)
    featureFile =  sessions{i,2};
    seizureFile = sessions{i,1};
    
    %get times used in training
    trainingTime = sort(cell2mat(cellfun(@(a) a(a(:,1)==i,2),sesID,'UniformOutput',false)'));
   % trainingTime = [];
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_rat(featureFile,seizureFile,rusTree,trainingTime,ops);
    outfil = [ClassifierFileOutputDir 'predict_' num2str(i) '.mat'];
    save(outfil,'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' outfil])
end

allUnlabel = [...;
    {'R:\DGregg\NeuralData\EDS_Cohort2\Induction Test\10-30-2023(14.35)\RHS_231030_143604\amplifier.lfp'}; ... ... ... ... ... ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHS1\11-17-2023(16.48)_BL\RHS_231117_164855\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHS1\11-17-2023(19.5)_BL\RHS_231117_190512\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-21-2023(15.1)\RHS_231121_150142\amplifier.lfp'         }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-21-2023(20.3)_BL\RHS_231121_200319\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-23-2023(16.18)_BL\RHS_231123_161828\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-23-2023(19.38)_BL\RHS_231123_193839\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-23-2023(2.2)_BL\RHS_231123_020246\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-23-2023(22.58)_BL\RHS_231123_225823\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(16.29)_BL\RHS_231124_162905\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(19.47)_BL\RHS_231124_194722\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(2.17)_BL\RHS_231124_021801\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(23.4)_BL\RHS_231124_230429\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-25-2023(16.17)_BL\RHS_231125_161755\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-25-2023(19.35)_BL\RHS_231125_193530\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-25-2023(2.21)_BL\RHS_231125_022136\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-25-2023(22.52)_BL\RHS_231125_225223\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-27-2023(16.28)_BL\RHS_231127_162847\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-27-2023(19.52)_BL\RHS_231127_195249\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-27-2023(23.16)_BL\RHS_231127_231611\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-28-2023(2.38)_BL\RHS_231128_023854\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(16.17)_BL\RHS_231202_161744\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-2-2023(19.36)_BL\RHS_231202_193639\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-3-2023(12.58)\RHS_231203_130000\amplifier.lfp'         }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-3-2023(16.18)_BL\RHS_231203_161833\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-3-2023(19.36)_BL\RHS_231203_193659\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-4-2023(16.19)_BL\RHS_231204_161947\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-4-2023(19.39)_BL\RHS_231204_193913\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-5-2023(13.2)\RHS_231205_130313\amplifier.lfp'          }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-5-2023(16.23)_BL\RHS_231205_162325\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-5-2023(19.43)_BL\RHS_231205_194335\amplifier.lfp'      }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-7-2023(13.4)\RHS_231207_130531\amplifier.lfp'          }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP2\12-8-2023(13.14)\RHS_231208_131450\amplifier.lfp'         }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP3\12-15-2023(13.8)\RHS_231215_130906\amplifier.lfp'         }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP3\12-16-2023(12.58)\RHS_231216_130000\amplifier.lfp'        }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP3\12-16-2023(16.23)_BL\RHS_231216_162338\amplifier.lfp'     }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-10-2023(12.58)\RHS_231010_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-11-2023(12.58)\RHS_231011_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-11-2023(18.43)\RHS_231011_184500\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-12-2023(13.15)\RHS_231012_131538\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-14-2023(12.58)\RHS_231014_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-16-2023(12.58)\RHS_231016_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-17-2023(12.58)\RHS_231017_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-18-2023(13.9)\RHS_231018_130954\amplifier.lfp'        }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-19-2023(13.5)\RHS_231019_130532\amplifier.lfp'        }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-20-2023(12.58)\RHS_231020_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-23-2023(12.58)\RHS_231023_130000\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\October\10-30-2023(18.53)\RHS_231030_185500\amplifier.lfp'       }; ...
    {'R:\DGregg\NeuralData\EDS_Cohort2\September\9-18-2023(18.36)\RHS_230918_183639\amplifier.lfp'      }; ...
    ]



for ii = 1:length(allUnlabel)
    
    s = dir(allUnlabel{ii});
dur = s.bytes/ops.nCh_featureFile/ops.Fs/2;



ts = 0: (dur-ops.durFeat);

%get prediction
estimateLabel =[];
dat1 =[];
for i = ts
    
    
    
    
    tim = i;
    features = ops.features(allUnlabel{ii},tim,ops);
    
    
    dat1 = [dat1;features];
    
    if mod(i,100)== 0
        [outpred,conf] = predict(rusTree,dat1);
        estimateLabel = [estimateLabel;outpred conf];
        dat1 =[];
    elseif i > (dur- mod(dur,100))
        [outpred,conf] = predict(rusTree,dat1);
        estimateLabel = [estimateLabel;outpred conf];
        dat1 =[];
    end
    
    
    
end
  outfil = [ClassifierFileOutputDir 'predict_u' num2str(ii) '.mat'];
    save(outfil,'estimateLabel')
end

%%
v=load(['R:\Analysis\SeizureForecasting\IHKA_rat_RF\Classification\predict_' num2str(19) '.mat'])
for i = 1:6
post(i,:) = nanconvn(v.estimateLabel(:,1)==i,k);
end

bar([nanmean(post(1:4,1:3600),2) nanmean(post(1:4,3601:7200),2) nanmean(post(1:4,7201:end),2)]')
%%

% load all predictions
topDir = 'R:\IHKA_Scharfman\prediction';
%topDir = ClassifierFileOutputDir;
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'predict') & ~contains(fils,'before');
fils  = fils(kp);

estimateLabel =[];trueLabel=[];time2seizure =[];confLabel =[];ID=[];
for i = 1:length(fils)
    
    v= load(fils{i});
    
    minLen = min(length(v.trueLabel),length(v.estimateLabel));
    trueLabel = [trueLabel;v.trueLabel(1:minLen)'];
    estimateLabel = [estimateLabel;v.estimateLabel(1:minLen,1)];
       confLabel = [confLabel;v.estimateLabel(1:minLen,2:end)];
    time2seizure = [time2seizure;v.time2seizure(1:minLen)'];
    ID = [ID;i*ones(minLen,1)];
    i
end
%%
clear AUC

AUCb = nan(max(ID),6,100);
for i = 1:max(ID)
     kp = ID ==i &trueLabel>0 ;
for j = 1:6
   if any(trueLabel(kp)==j)
[X,Y,T,AUC(i,j)] = perfcurve(trueLabel(kp)==j,confLabel(kp,j),1); 
sen(i,j) = nanmean(estimateLabel(trueLabel==j&kp)==j);
sp(i,j) = nanmean(estimateLabel(trueLabel~=j&kp)~=j);

%get null

confLabelR = confLabel(kp,j);
%for b = 1:100
  %  confLabelR = confLabelR(randsample(1:length(confLabelR),length(confLabelR)));
 %   [X,Y,T,AUCb(i,j,b)] = perfcurve(trueLabel(kp)==j,confLabelR,1); 

%end
    
   else
       AUC(i,j) = nan;
       
   end

end
end
%%
kp = AUC(:,4)>.5;
figure
plot(AUC(kp,:)','color',[.7 .7 .7])
hold on
plotMeanSEM(1:6,AUC(kp,:),'k')
plotMeanSEM(1:6,nanmean(AUCb(kp,:,:),3),'r')
%%
C = confusionmat(trueLabel(trueLabel>0),estimateLabel(trueLabel>0));

figure
imagesc(C./sum(C,2),[0 1])
set(gca,'yticklabel',{'3hrs','1hr','10min','10s','sz','post'})
set(gca,'xticklabel',{'3hrs','1hr','10min','10s','sz','post'})
xlabel('Predicted class')
ylabel('Real class')

set(gca,'ydir','normal')

%%
close all
figure
k  = gaussian2Dfilter([1000 1],0);
for i = 1:5
    d = double(estimateLabel==i);
    %d = nanconvn(d,k);
semilogx(fliplr(avghist(time2seizure,d,-43000:0)))
hold on
end

%%

figure
k  = gaussian2Dfilter([1000 1],3);
for i = 1:5
    d = double(estimateLabel==i);
   
semilogx(fliplr(avghist(time2seizure,d,-86400:0,@numel)))
hold on
end


%%



close all
figure
k = gaussian2Dfilter([1 100],[ 1 50]);
%estimateLabel1 = estimateLabel;
%estimateLabel1(inTrainingSet) = nan;
for i = 1:6
    subplot(6,1,i)
   
    hold on
    plot([seizure_start seizure_start],[0 1],'--','color','r')
     plot(nanconvn(estimateLabel==i,k'),'k')
end



%%


%get ROC''
clear FP TP AUC T AUCtmp lb hb
ix=1;
for i=[1 2 3 4 5 6]
    tmp =  group(ops.N+1:end) ;
%tmp(tmp==2) = 3;
[FP{ix},TP{ix},T{ix},AUC(ix)] = perfcurve(tmp,score(:,i),i);


% get random

for j = 1:100
    tmpt = tmp(randsample(1:length(tmp),length(tmp)));
    [~,~,~,AUCt(j)] = perfcurve(tmpt,score(:,i),i);
end
AUCtmp(ix) = mean(AUCt);
lb(ix) = prctile(AUCt,1);
hb(ix) = prctile(AUCt,99);
ix = ix+1;
end

hb = hb-AUCtmp;
lb = lb-AUCtmp;


%%

    %%
    
     load('E:\data\IHKA\classifier1.mat')
     
    %%
   
   
  
  
    [pred,score] = predict(rusTree,training);

actual_Y = group;
predicted_Y = pred;

%%

%%
close all

figure
x=-([18000,3600 600 10 2 1]);
y= AUC;
semilogx(x,y,'o-')
hold on
xx = [x';flipud(x')];
zz = [AUCtmp'+lb';flipud(AUCtmp')+hb'];



patch(xx,zz,'k','facealpha',.25,'edgecolor','none')
semilogx(x,AUCtmp,'--','color','k')
xlim([-18000 0])
%%


%%

ok = cell2mat(dat');
X  = tsne(nanzscore(ok));
group = cellfun(@(a,b) a*ones(length(b),1),num2cell(1:6),dat,'UniformOutput',false);
group = cell2mat(group');
[~,b] = histc(group,1:6);
kp = ~any(isnan(ok),2);
%%
close all
k  = gaussian2Dfilter([100 100],5);
ax  = tight_subplot(1,6);
for  ii = 1:6
    axes(ax(ii))
bin = histcn([X(b(kp)==ii,1),X(b(kp)==ii,2)],-50:50,-50:50)/sum(b(kp)==ii);
PP(ii,:) = bin(:);
imagesc(nanconvn(bin,k))
end

%%

for i = 1:6
    
    for j = 1:6
        
     kl(i,j) =   KLDiv(PP(i,:),PP(j,:));
        
    end
end