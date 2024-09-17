function idx_r=fcn_get_rising_edge(y)
    yy=y;%yy(isnan(y))=0;yy(~isnan(y))=1;      
    dy=diff([0,yy]);
    idx_r=find(dy==1);
end

