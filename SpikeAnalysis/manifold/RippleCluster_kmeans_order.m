clear all
topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\AAC\Acute';
%topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\PV';
%topDir = 'R:\AAC\AAC_DataForSam\ArchSessions\CCK';
plotIt = false;
d = dir(topDir);
d = {d(cell2mat({d.isdir})).name};
dirs = d(3:end);

dirs = cellfun(@(a) [topDir filesep a],dirs,'uni',0);
kp = ~contains(dirs,'figures');
dirs = dirs(kp);
clear p_comp uRate_comp spont_order stim_order
%%
% loop over directories

for i = 1:length(dirs)

    cd(dirs{i})
    fils = getAllExtFiles(pwd,'mat',0);

    kp = contains(fils,'celltypes') |  contains(fils,'ripples') |  contains(fils,'optoStim')  |  contains(fils,'spikes') | contains(fils,'cell_metrics');
    fils = fils(kp);
    %load all data

    %load data
    clear allcelltypes ripples spikes optoStim cell_metrics
    for j = 1:length(fils)
        load(fils{j})
    end

    if exist('allcelltypes')
        ispyr =cellfun(@(a) isstr(a) && contains(a,'pyr'),allcelltypes);

    else
        ispyr =cellfun(@(a) isstr(a) && contains(a,'Pyr'),cell_metrics.putativeCellType);
    end
    in = InIntervals(ripples.peaks,optoStim.timestamps);

    [binnedPop_post,bin_times]=populationMatrix(spikes,0,.1,1,ripples.timestamps(:,1));
    binSpk = squeeze(binnedPop_post)';
    zpop = zscore(binSpk(:,ispyr));
   uRate = sum(zpop(:,:),2);

    zpop1 = zpop./(repmat(uRate-min(uRate(uRate<0)),1,size(zpop,2)));


    kp_noinf = ~any(isinf(zpop1),2);
   
    % ok = nan(size(zpop1,1),2);
    %  ok(kp,:) = tsne(zpop1(kp,:));
    rip_ts = ripples.timestamps;

    % do kmean with cross validation


    max_k = 50;
    avg_silhouette = zeros(max_k-1, 1);  % No silhouette for k=1

    idxB = find(in& kp_noinf);   % trials to test
    idxA = setdiff([1:length(in)]',[idxB;find(~kp_noinf)]);   % manifold-defining trials
    rIDXA = idxA(randsample(length(idxA),length(idxA)));
    nBlock = 10;
    nShuff = 10;  % number of shuffles for controls
    blidx = 1:floor(length(rIDXA)/nBlock):length(rIDXA);
    [n,bl] = histc(1:length(rIDXA),blidx);

    %define k-fold
    for b = 1:nBlock

        manA = rIDXA(bl ~=b);
        testA = rIDXA(bl ==b);
        nB = length(idxB);
        varN = num2cell(num2str([1:size(zpop1,2)]'),2);
        [old_embeddinga,embedding]  = run_umap(zpop1(manA,:), 'save_template_file', 'TEST.MAT','verbose', 'none','plot', 'none','see_training'  ,false, 'parameter_names',varN,'min_dist',.05,'n_neighbors',10);
        [new_embeddinga] = run_umap(zpop1(testA,:), 'template_file', 'TEST.MAT','verbose', 'none','plot', 'none','see_training'  ,false,'parameter_names', varN);
        [new_embeddingb] = run_umap(zpop1(idxB,:), 'template_file', 'TEST.MAT','verbose', 'none','plot', 'none','see_training'  ,false,'parameter_names', varN);
        %umap_embedded = [old_embeddinga;new_embeddinga;new_embeddingb];


        umap_embedded = nan(size(zpop1,1),2);
        umap_embedded(manA,:) = old_embeddinga;
        umap_embedded(testA,:) = new_embeddinga;
        umap_embedded(idxB,:) = new_embeddingb;



        for k = 5:max_k
            % Cluster using k-means
            idx = kmeans(umap_embedded, k, 'Replicates', 5, 'Display', 'off');

            % Compute silhouette values
            s = silhouette(umap_embedded, idx);

            % Average silhouette score
            avg_silhouette(k-4) = nanmean(s);
        end


        [~, best_k] = max(avg_silhouette);
        optimal_k = best_k + 4;
 
        ix = kmeans(umap_embedded,optimal_k);

        p_comp{i,b} = [];uRate_comp{i,b} =[];
        for ii = 1:optimal_k
            % define probabliity of new data
            kp_base = ismember(1:size(umap_embedded,1),manA)' & ix==ii;
            kp_held_out = ismember(1:size(umap_embedded,1),testA)' & ix==ii;
            kp_stim =  ismember(1:size(umap_embedded,1),idxB)' & ix==ii;
          
            

            if sum(kp_base)>5
                mu =  nanmean(  umap_embedded(kp_base,:),1);
                sigma =  nancov(  umap_embedded(kp_base,:))+1e-10;
                new_data = umap_embedded(kp_stim,:);
                old_data = umap_embedded(kp_held_out,:);
                p_comp{i,b}(ii,1) = mean(log(mvnpdf(new_data, mu, sigma)));
                p_comp{i,b}(ii,2) = mean(log(mvnpdf(old_data, mu, sigma)));
            else

                p_comp{i,b}(ii,:) = [nan nan];
            end
        end




        uRate_comp{i,b}(:,1)  =accumarray(ix(kp_noinf&in),uRate(kp_noinf&in),[optimal_k 1],@nanmean,nan);
        uRate_comp{i,b}(:,2) =accumarray(ix(kp_noinf&~in),uRate(kp_noinf&~in),[optimal_k 1],@nanmean,nan);



        nBin = 5;
        [binnedPop_post,bin_times]=populationMatrix(spikes,0,.1,nBin,ripples.timestamps(:,1),'zscore',true);
        binnedPop_post = binnedPop_post(ispyr,:,:);
        spont_order{i,b} = []; stim_order{i,b} = [];
        for ii = 1:optimal_k
            %figure
            %ax  = tight_subplot(2,1);
            %axes(ax(1))
            kp_base = ismember(1:size(umap_embedded,1),manA)' & ix==ii;
            kp_held_out = ismember(1:size(umap_embedded,1),testA)' & ix==ii;
            kp_stim =  ismember(1:size(umap_embedded,1),idxB)' & ix==ii;
            
            if sum(kp_base)>5
                %get sequence order
                [a,b_order ] =max(nanmean(binnedPop_post(:,:,kp_base),3),[],2);
                in_clust = a > .15;

                %get held out rate
                tmp = nanmean(binnedPop_post(:,:,kp_held_out),3);

                rip_clust_spont{i,b}{ii} = tmp;
                tmp1 = nan(nBin,nBin);
                for oo = 1:nBin
                    tmp1(:,oo) = accumarray(b_order(in_clust),tmp(in_clust,oo),[nBin 1],@nanmean,nan);

                end

                spont_order{i,b}(:,:,ii) = tmp1;

                %  imagesc(tmp(,[-1 1]*.5)

                % axes(ax(2))


                tmp = nanmean(binnedPop_post(:,:,kp_stim),3);

                rip_clust_stim{i,b}{ii} = tmp;
                tmp1 = nan(nBin,nBin);
                for oo = 1:nBin
                    tmp1(:,oo) = accumarray(b_order(in_clust),tmp(in_clust,oo),[nBin 1],@nanmean,nan);

                end

                stim_order{i,b}(:,:,ii) = tmp1;

                %  imagesc(tmp(kp1,:),[-1 1]*.5)
                %  set(gcf,'position',[  1          41        1920        1083])
                %  waitforbuttonpress
                %  close all
            else
                stim_order{i,b}(:,:,ii) = nan(nBin,nBin);
                spont_order{i,b}(:,:,ii) = nan(nBin,nBin);
                rip_clust_stim{i,b}{ii} = [];
                rip_clust_spont{i,b}{ii} = [];

            end
        end
    end
    %
    % if plotIt
    %     close all
    %
    %
    %     figure
    %
    %     col = linspecer(optimal_k,'jet');
    %     for ii = 1:optimal_k
    %
    %         plot(ok(ix==ii& ~in,1),ok(ix==ii&~in,2),'.','color',col{ii})
    %
    %         hold on
    %     end
    %     plot(ok(in,1),ok(in,2),'x','color','k')
    %     plot(ok(in,1),ok(in,2),'o','color','k')
    %
    %     recenter = ok;
    %
    %     for ii = 1:optimal_k
    %
    %         recenter(ix==ii,:)  =   recenter(ix==ii,:) - nanmean(  recenter(ix==ii & ~in,:));
    %     end
    %     figure
    %     plot(recenter(~in&kp,1),recenter(~in&kp,2),'.','color',[.7 .7 .7])
    %     hold on
    %     plot(recenter(in&kp,1),recenter(in&kp,2),'x','color','k')
    %     plot(recenter(in,1),recenter(in,2),'o','color','k')
    %
    %
    %
    %
    %
    %     figure
    %     plot(0:20,0:20)
    %     hold on
    %     plot( uRate_comp{i}(:,2),uRate_comp{i}(:,1),'.')
    %     xlabel('no stim')
    %     ylabel('stim')
    % end
    %
end
%%
p_comp1 = [];
for i = 1:size(p_comp,1)
    tmp = p_comp(i,:);
    p_comp1  = [p_comp1 ; cell2mat(cellfun(@(a) a(:,1),tmp,'UniformOutput',false)') ...
        cell2mat(cellfun(@(a) a(:,2),tmp,'UniformOutput',false)')];
end



uRate_comp1 = [];
for i = 1:size(p_comp,1)
    tmp = uRate_comp(i,:);
    uRate_comp1  = [uRate_comp1 ; cell2mat(cellfun(@(a) a(:,1),tmp,'UniformOutput',false)') ...
        cell2mat(cellfun(@(a) a(:,2),tmp,'UniformOutput',false)')];
end



%%

figure
plot(p_comp1(:,2),p_comp1(:,1),'o')
hold on
xlabel('LL spont')
ylabel('LL stim')
plot(-7:.5:-.5,-7:.5:-.5)
axis square
%%

figure
plot(uRate_comp1(:,2),uRate_comp1(:,1),'o')
hold on
xlabel('rate mod spont')
ylabel('rate mod stim')
plot(-20:.5:35,-20:.5:35)
axis square
%%
stim_ok = cellArrayTo3D(cellfun(@(a) nanmean(a,3),stim_order(:,:),'UniformOutput',false));
spont_ok = cellArrayTo3D(cellfun(@(a) nanmean(a,3),spont_order(:,:),'UniformOutput',false));




i


close all
    imagesc(nanmedian( spont_ok,3))
title('Spont')
figure
imagesc(nanmedian(stim_ok,3))
shg

title('Stim')


%%



[p,h] = signtest(uRate_comp1(:,1),uRate_comp1(:,2))

[p,h] = signtest(p_comp1(:,1),p_comp1(:,2))

%%

stim_order1 = cellArrayTo3D(cellfun(@(a) nanmean(a,3),stim_order,'uni',0));
spont_order1 = cellArrayTo3D(cellfun(@(a) nanmean(a,3),spont_order,'uni',0));
%%
figure
imagesc(nanmean(spont_order1,3))
figure
imagesc(nanmean(stim_order1,3))

%%
for k = 1:10
    order_rip_spont{k} =[];
    order_rip_stim{k} =[];

    for i = 1:length(rip_clust_spont)
        for j = 1:length(rip_clust_spont{i})
            if size( rip_clust_spont{i}{j},1) == size( rip_clust_stim{i}{j},1)
                tmp1 = rip_clust_spont{i}{j};
                tmp2 = rip_clust_stim{i}{j};

                [~,b] = max(tmp1,[],2);




                order_rip_spont{k} = [order_rip_spont{k} ; tmp1(b==k,:)];
                order_rip_stim{k} = [order_rip_stim{k} ; tmp2(b==k,:)];

            end

        end
    end
end


%%
close all

ax  = tight_subplot(10,1);
ii=1;
for i = 1:1:10
    %    axes(ax(ii))
    rip_mat_spont(i,:) = (nanmean(cell2mat(order_rip_spont(i)')));
    rip_mat_stim(i,:) =(nanmean(cell2mat(order_rip_stim(i)')));
    %   ii=ii+1;
end

%%
close all
imagesc(rip_mat_spont)
figure
imagesc(rip_mat_stim)

