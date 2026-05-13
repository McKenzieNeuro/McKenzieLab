
function aligned = iterative_procrustes(raw_cell, nIter)
    n=length(raw_cell);
    aligned=raw_cell;
    for k=1:n, aligned{k}=aligned{k}-mean(aligned{k},1); end
    for iter=1:nIter
        mat=cat(3,aligned{:}); mu=mean(mat,3); mu=mu-mean(mu,1);
        for k=1:n
            [~,aligned{k}]=procrustes(mu,aligned{k},...
                'Scaling',false,'Reflection',false);
        end
    end
end