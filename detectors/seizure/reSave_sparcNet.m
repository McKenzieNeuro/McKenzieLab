function reSave_sparcNet(fil)


Fs=200;ww=2*Fs;


load(fil)
[xx,num_points]=size(data);
num_seg=ceil(num_points/ww);

featue_file=strrep(fil,'.mat','_score.csv');
outfil = strrep(featue_file,'csv','mat');

Y=cell2mat(table2cell(readtable(featue_file)));
Y=[repmat(Y(1,:),2,1);Y];
Y=[Y;repmat(Y(end,:),num_seg-size(Y,1),1)];
[~,y0]=max(Y,[],2);y0=y0-1;


% parameters
Fs=200;ww=2*Fs;
params.movingwin=[4,2];params.tapers=[2,3];params.fpass=[.5,20];params.Fs=Fs;
thr_cpd=.1;sgolay_para=10;num_words=20;num_clusters=30;thr_dur=5;
tmp=load('bsd_lut.mat');
K=10;lut_boost=tmp.lut(1:K,:);

num_seg=ceil(size(data,2)/ww);
% spectrogram
[sdata,~,sfreqs]=fcn_computeSpec(data,params);
for j=1:size(sdata,1)
    s=sdata{j,2};
    s=[s(:,1),s(:,1:end-1),repmat(s(:,end-1),1,num_seg-size(s,2))];
    sdata{j,2}=s;
end
stimes=(0:2:(num_seg-1)*2)+1;

% cpd
s=cell2mat(sdata(:,2));
P=mean(pow2db(s+eps));P=smooth(P,sgolay_para,'sgolay')';P(P>25)=25;P(P<-10)=-10;pow=P;
idx_cp=findchangepts(P,'Statistic','mean','MinThreshold',thr_cpd*var(P));
idx_cp=unique([1,idx_cp,length(P)]);
idx_cp1=idx_cp(1:end-1);
idx_cp2=idx_cp(2:end)-1;idx_cp2(end)=idx_cp(end);
idx_cpc=floor((idx_cp1+idx_cp2)/2);

% bow cluster
X=pow2db(cell2mat(sdata(:,2))+eps)';
rng('default');
y=kmeans(X,num_words);
num_seg=length(idx_cp)-1;
bow=NaN(num_seg,num_words);
for k=1:num_seg-1
    z=y(idx_cp(k):idx_cp(k+1)-1);
    bow(k,:)=hist(z,1:num_words);
end
z=y(idx_cp(num_seg):idx_cp(num_seg+1));
bow(num_seg,:)=hist(z,1:num_words);
if num_seg<num_clusters
    num_clusters = num_seg;
end
bow_n=bow./repmat(sum(bow,2),1,size(bow,2));
rng('default');
[idx_clusters,medoids]=kmedoids(bow_n,num_clusters,'Distance',@fcn_distChiSq);
idx_medoids=NaN(num_clusters,1);
idx_members=NaN(num_seg,1);
for k=1:num_clusters
    kk=find(ismember(bow_n,medoids(k,:),'rows'));
    idx_medoids(k)=kk(1);
    idx_members(idx_clusters==k)=kk(1);
end
idx_medoids=sort(idx_medoids);

% post-sparcnet processing
y1=NaN(num_seg,1);y2=y1;
for k=1:length(idx_cpc)
    idx_seg=idx_cpc(k);
    [~,yh]=max(Y(idx_seg,:));
    pred_org=yh-1;
    
    pred=pred_org;
    if pred>0
        % parse eeg
        tc=(((idx_seg-1)*2+1)*Fs)+1;
        seg=fcn_parse_eeg(data,tc,10*Fs);
        
        % compute features
        fff=fcn_compute_bgr10(seg,Fs);
        
        % eliminate "other" via boosting
        is_other=0;
        for kk=1:size(lut_boost,1)
            zzz=sign(lut_boost{kk,3})*fff(kk);
            thr=lut_boost{kk,4};
            if zzz<thr
                is_other=1;
                break;
            end
        end
        if is_other==1
            pred=0;
        end
    end
    
    % part #1: smooth in cp segments
    y1(idx_cp1(k):idx_cp2(k))=pred_org;
    
    % part #2: eliminate other via boosting
    y2(idx_cp1(k):idx_cp2(k))=pred;
end

% part #3: check amplitude / power
y3=y2;
y3(pow<=-10|pow>=20)=0;

% part #4: check event duration
y4=y3;
for k=2:6
    yy=(y4==(k-1))';
    idx_r=fcn_get_rising_edge(yy);
    idx_f=fcn_get_falling_edge(yy);
    dur=idx_f-idx_r+1;
    for kk=1:length(dur)
        if dur(kk)<thr_dur
            try
                y4(idx_r(kk):idx_f(kk))=y4(idx_r(kk)-1);
            catch err
                y4(idx_r(kk):idx_f(kk))=y4(idx_f(kk)+1);
            end
        end
    end
end

% part #5: smooth in bow clusters
y5=y4;
for k=1:num_clusters
    vv=hist(y5(idx_cpc(idx_members==idx_medoids(k))),0:5);
    [~,yy]=max(vv);
    y5(idx_cpc(idx_members==idx_medoids(k)))=(yy-1);
end
for k=1:length(idx_cpc)
    y5(idx_cp1(k):idx_cp2(k))=y5(idx_cpc(k));
end
y5(pow<=-10|pow>=20)=0;
for k=2:6
    yy=(y5==(k-1))';
    idx_r=fcn_get_rising_edge(yy);
    idx_f=fcn_get_falling_edge(yy);
    dur=idx_f-idx_r+1;
    for kk=1:length(dur)
        if dur(kk)<thr_dur
            try
                y5(idx_r(kk):idx_f(kk))=y5(idx_r(kk)-1);
            catch err
                y5(idx_r(kk):idx_f(kk))=y5(idx_f(kk)+1);
            end
        end
    end
end
%%
% export
pred=y0;pred_plus=y5;
save(outfil,'pred','pred_plus','Y')

end

