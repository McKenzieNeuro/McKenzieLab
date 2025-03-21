function ff=fcn_compute_bgr10(seg,Fs)
    params.movingwin=[4,.5];params.tapers=[2,3];params.fpass=[.5,20];params.Fs=Fs;       

    nn=size(seg,2);
    gg=reshape(1:16,4,4)'; % LL RL LP RP
    gl={'LL','RL','LP','RP'};
    
    [Sdata,~,sfreqs]=fcn_computeSpec(seg,params);
  
    dur=nn/Fs;
    stimes=0:.5:dur;
       
    % f1:'nleo_std-LP-2sec'
    w=2;r='LP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    seg_=seg(gg(ismember(gl,r),:),tt1:tt2);
    ff=NaN(1,size(seg_,1));
     
    for ii=1:size(seg_,1)
        x=seg_(ii,:);
        x_nleo=general_nleo(x);
        ff(ii)=std(x_nleo);
    end
    f1=mean(ff);

    % f2:line_length-RP-2sec
    w=2;r='RP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    seg_=seg(gg(ismember(gl,r),:),tt1:tt2);
    ff=NaN(1,size(seg_,1));
    for ii=1:size(seg_,1)
        x=seg_(ii,:);
        ff(ii)=mean(abs(diff(x)));
    end
    f2=mean(ff);

    % f3:'power_ratio_theta2alpha-RP-10sec-5%'
    w=10;r='RP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    s_tt1=find(round(stimes*10)/10==round(10*tt1/Fs)/10);
    s_tt2=find(round(stimes*10)/10==round(10*tt2/Fs)/10);
    spec=Sdata{ismember(gl,r),2};
    spec=[repmat(spec(:,1),1,4),spec(:,1:end-1),repmat(spec(:,end-1),1,5)]; 
     
    ss=spec(:,s_tt1:s_tt2);

    pp_theta=bandpower(ss,sfreqs,[4 8],'psd');
    pp_alpha=bandpower(ss,sfreqs,[8 12],'psd');
   
    pr_t2a=pp_theta./(pp_alpha+eps);
     
    
    f3=prctile(pr_t2a,5);

    %  f4:'shannon_entropy-LL-2sec'
    w=2;r='LL';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    
    seg_=seg(gg(ismember(gl,r),:),tt1:tt2);
    ff=NaN(1,size(seg_,1));
     
    for ii=1:size(seg_,1)
        x=seg_(ii,:);
        ff(ii)=mean(pentropy(x,Fs));
    end
    
    f4=mean(ff);
   
    % f5:power_beta-LP-2sec
    w=2;r='LP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    s_tt1=find(round(stimes*10)/10==round(10*tt1/Fs)/10);
    s_tt2=find(round(stimes*10)/10==round(10*tt2/Fs)/10);
    spec=Sdata{ismember(gl,r),2};
    spec=[repmat(spec(:,1),1,4),spec(:,1:end-1),repmat(spec(:,end-1),1,5)];        
    ss=spec(:,s_tt1:s_tt2);
 
    pp_beta=bandpower(ss,sfreqs,[12 18],'psd');
    f5=mean(pp_beta);
  
    %  f6:'power_ratio_theta2alpha-LP-10sec-95%'
    w=10;r='LP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    s_tt1=find(round(stimes*10)/10==round(10*tt1/Fs)/10);
    s_tt2=find(round(stimes*10)/10==round(10*tt2/Fs)/10);
    spec=Sdata{ismember(gl,r),2};
    spec=[repmat(spec(:,1),1,4),spec(:,1:end-1),repmat(spec(:,end-1),1,5)];        
    ss=spec(:,s_tt1:s_tt2);

    pp_theta=bandpower(ss,sfreqs,[4 8],'psd');
    pp_alpha=bandpower(ss,sfreqs,[8 12],'psd');
    pr_t2a=pp_theta./(pp_alpha+eps);
    f6=prctile(pr_t2a,95);

    %  f7:'relative_power_alpha-LP-10sec-95%'
    w=10;r='LP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    s_tt1=find(round(stimes*10)/10==round(10*tt1/Fs)/10);
    s_tt2=find(round(stimes*10)/10==round(10*tt2/Fs)/10);
    spec=Sdata{ismember(gl,r),2};
    spec=[repmat(spec(:,1),1,4),spec(:,1:end-1),repmat(spec(:,end-1),1,5)];        
    ss=spec(:,s_tt1:s_tt2);

    pp_alpha=bandpower(ss,sfreqs,[8 12],'psd');
    pp_total=bandpower(ss,sfreqs,[sfreqs(1) sfreqs(end)],'psd');
    
    rp_alpha=pp_alpha./(eps+pp_total);
    f7=prctile(rp_alpha,95);

    % f8:'power_ratio_delta2alpha-RL-10sec-50%'
    w=10;r='RL';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    s_tt1=find(round(stimes*10)/10==round(10*tt1/Fs)/10);
    s_tt2=find(round(stimes*10)/10==round(10*tt2/Fs)/10);
    spec=Sdata{ismember(gl,r),2};
    spec=[repmat(spec(:,1),1,4),spec(:,1:end-1),repmat(spec(:,end-1),1,5)];        
    ss=spec(:,s_tt1:s_tt2);

    pp_delta=bandpower(ss,sfreqs,[1 4],'psd');
    pp_alpha=bandpower(ss,sfreqs,[8 12],'psd');
    pr_d2a=pp_delta./(pp_alpha+eps);
    f8=prctile(pr_d2a,50);

    % f9:'nleo_std-LP-6sec'
    w=6;r='LP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    seg_=seg(gg(ismember(gl,r),:),tt1:tt2);
    ff=NaN(1,size(seg_,1));
    for ii=1:size(seg_,1)
        x=seg_(ii,:);
        x_nleo=general_nleo(x);
        ff(ii)=std(x_nleo);
    end
    f9=mean(ff);

    % f10:'nleo_std-RP-6sec'
    w=6;r='RP';
    tt1=nn/2-w*Fs/2+1;tt2=tt1+w*Fs-1;
    seg_=seg(gg(ismember(gl,r),:),tt1:tt2);
    ff=NaN(1,size(seg_,1));
    for ii=1:size(seg_,1)
        x=seg_(ii,:);
        x_nleo=general_nleo(x);
        ff(ii)=std(x_nleo);
    end
    f10=mean(ff);


    ff=[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10];
end