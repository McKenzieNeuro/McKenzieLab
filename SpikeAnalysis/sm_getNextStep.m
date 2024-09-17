function pred = sm_getNextStep(mat,back)

pred = nan(size(mat));
for i = back+1:size(mat,2)
    
    a=ltv_adjacency(mat(:,i-back:i-1)');
    pred(:,i) = a*mat(:,i-1);

end