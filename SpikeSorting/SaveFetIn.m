
function SaveFetIn(FileName, Fet, BufSize)

if nargin<3 | isempty(BufSize)
    BufSize = inf;
end

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
    formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');
fprintf(outputfile, '%d\n', nFeatures);

if isinf(BufSize)
    
    temp = [round(100* Fet(:,1:end-1)) round(Fet(:,end))];
    fprintf(outputfile,formatstring,temp');
else
    nBuf = floor(size(Fet,1)/BufSize)+1;
    
    for i=1:nBuf
        BufInd = [(i-1)*nBuf+1:min(i*nBuf,size(Fet,1))];
        temp = [round(100* Fet(BufInd,1:end-1)) round(Fet(BufInd,end))];
        fprintf(outputfile,formatstring,temp');
    end
end
fclose(outputfile);
end