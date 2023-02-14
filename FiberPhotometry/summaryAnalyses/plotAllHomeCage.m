masterDir = 'R:\McKenzieLab\DANEHippocampalResponse';

masterDir = 'Y:\DANEHippocampalResponse';
fils = getAllExtFiles(masterDir,'mat',1);

kp = cellfun(@any,regexp(fils,'Novel Env'));

fils = fils(kp);

[dirs] = cellfun(@fileparts,fils,'uni',0);

dirs =  unique(dirs);
%%
k  = gaussian2Dfilter([100 100],[1 1]);
ixx = 1;
for i = 30:length(dirs)
    i
    cd(dirs{i})
    

    load('contextTransition1.mat');

    if any(cellfun(@any,regexp(data(:,1),'WB')))
        if exist('sessiondata.mat')

            %load tracking data
            load('sessiondata.mat');

            %load context edges
            %   load('contextTransition1.mat');

            kp = cellfun(@any,regexp(sessiondata.behavior.position.context,'WB'));
            LE =  sessiondata.behavior.position.left_ear_cor;
            LE = LE(kp,:);
            ts = sessiondata.behavior.ts_video;
            ts_neural  = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
            dsNeural = interp1(ts_neural,sessiondata.neural.signal_DFoF,ts);


            ts = ts(kp);
            dsNeural = double(dsNeural(kp));

            bins1 = -1:.1:1.2;
            bins2 = -1:.1:1.2;
            [occ,~,~,ix] = histcn(LE,bins1,bins2);
            inbin = all(ix>0,2);
            tmp = accumarray(ix(inbin,:),dsNeural(inbin),[length(bins1),length(bins2)],@nanmean,nan);
            tmp(occ<10) = nan;
            ratemap(:,:,ixx) = nanconvn(tmp,k,'nanout',true);
            sess{ixx} = dirs{i};
            ixx = ixx+1;
        end

    end
end

%%
n = 100;
x0 = 6;
y0 = 15;

for i = 1:100
R = i;

theta=linspace(0,360,n+1); theta(end)=[];
x=R*cosd(theta)+x0;
y=R*sind(theta)+y0;

[~,~,~,ix] = histcn([x' y'],1:length(bins1),1:length(bins2));
ix = unique(ix,'rows');
ix = ix(all(ix>0,2),:);
ix_1 = repmat(ix(:,2),size(ratemap,3),1);

ix_2 = repmat(ix(:,1),size(ratemap,3),1);

ix_3 = repmat([1:size(ratemap,3)]',size(ix,1),1);

ok  =sub2ind(size(ratemap),ix_1,ix_2,ix_3);
sig_Rc{i} = (ratemap(ok));
end
