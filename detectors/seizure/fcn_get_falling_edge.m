function idx_f=fcn_get_falling_edge(y)
    yy=y;%yy(isnan(y))=0;yy(~isnan(y))=1;
    dy=diff([yy,0]);
    idx_f=find(dy==-1);
end