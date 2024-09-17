function data=fcn_parse_eeg(data,tc,w)%(eeg_file,tc,Fs,w)

    [num_ch,num_pt]=size(data); 
    t1=max(1,tc-w/2+1);t2=min(num_pt,tc+w/2);
    seg=data(:,t1:t2);

    data=zeros(num_ch,w);
    L2=num_pt;
    i1=tc-w/2+1;i2=tc+w/2;

    % reasoning
    if i1<1&&i2<L2 % ---xxx 
        %disp(' Early')
        bb=w;
        aa=bb-size(seg,2)+1;
        data(:,aa:bb)=seg;
    elseif i1<1&&i2>=L2 % --xx--
        %disp(' Short')
        aa=1-i1+1;
        bb=aa+L2-1;
        data(:,aa:bb)=seg;

    elseif i1>=1 && i2<L2  % xxxxxx
        %disp(' Regular')
        data=seg;

    elseif i1>=1 && i2>=L2 % xxx---
        %disp(' Late')
        aa=1;
        bb=size(seg,2);
        data(:,aa:bb)=seg;
    else
        keyboard;
    end
end
