
function [prob,est] = decodeSimple(spk_times,linearPos)
 
%spike times is a cell array of each cells spike times

%linearPos column 1 =ts, column 2 = position at each time stamp

%define temporal and spatial bins
binsize = .2; %200ms

allpos = linearPos(:,2);
spatbin = 10; %pixel resolution
pos = min(allpos):spatbin:max(allpos); %bins for ratemapt in pixels



%sampling rate
fs = mode(diff(linearPos(:,1)));

%define time stamps for training and decoding
%this ensures training on half the data and decoding on the other half
ts=[min(linearPos(:,1)):binsize:max(linearPos(:,1))];
ev=(fs/binsize)-mod(length(ts),(fs/binsize));
ts=[ts (ts(end)+binsize):binsize:(ts(end)+ev*binsize)];
temp=reshape(ts,(fs/binsize),[])';
train=temp(1:2:end,:);
train=any(histc(train(:),ts),2);
decode=temp(2:2:end,:);
decode=any(histc(decode(:),ts),2);


%define position at each time
if binsize <= fs^-1
    
    %interpret position at higher sampling rate
    linearPos = interp1(linearPos(:,1),linearPos(:,2),ts);
   
else
    
    %downsample position to average within binsize
    linearPos = avghist(linearPos(:,1),linearPos(:,2),ts);
    
end

%bin space
[occ,bin]=histc(linearPos,pos);
occ = occ*binsize;


%set up smoothing kernel
kernel = pdf('Normal', -slide:slide+1, 0, 2);
%kernel=1;

%loop through cells for getting training posterior distribution (ratemaps)
for j=1:size(cel,1)
    
    
    spkstemp=histc(spk_times{j},ts);
  
    
    
    %delete spikes where position is not defined
    spkstemp(bin==0)=nan;
    
   
      
      
    spks(:,j)=spkstemp;
    
    temp=avghist(linearPos(train),spkstemp(train),pos)./occ';
    temp=temp(1:end-1);
    
    %convolve rate map
    recField(:,j)=nanCircularConvn(temp, kernel); % great function! I am sending it your way
    
   
end
    
%calculate posterior probs
for b=1:spatbin
    
    prob(:,b)=nanprod(((repmat(binsize*recField(b,:),size(spks,1),1).^spks)).*exp(repmat(-binsize*recField(b,:),size(spks,1),1)),2);
  
    
  
end

%normalize
prob=prob./repmat(sum(prob,2),1,size(prob,2));

%nan out the training time stamps
prob(train,:)=nan;

%get position estimate

[~,est] = max(prob,[],2);

%delete conditions with no spikes, with uniform posterior and during training
est(nanvar(prob,[],2)==0) = nan;
est(all(isnan(prob),2)) = nan;
est(all(spks==0,2)) = nan;
end


function [hist_y, n_y,hist_ySEM]= avghist(x,y,edges)
%returns a histogram of the average values of y based upon binning the data
%in x
if isempty(x)
    hist_y=[];
    n_y=[];
    return;
else
    
    edges=edges(:);
    if length(y)<length(x)
        x=x(1:length(y));
        
    end
    
    [n,bin]=histc(x,edges);
    
    if all(n==0)
        hist_y=[];
        n_y=[];
        return;
    end
    bin=bin(:);
    y=y(:);
    
    hist_y=accumarray(bin(bin~=0),y(bin~=0),[length(edges) 1],@nanmean,nan)';
    hist_ySEM=accumarray(bin(bin~=0),y(bin~=0),[length(edges) 1],@SEM,nan)';
    n_y=accumarray(bin(bin~=0),y(bin~=0),[length(edges) 1],@(a) nnz(~isnan(a)),nan)';
    
end
end