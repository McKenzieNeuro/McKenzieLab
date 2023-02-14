function addOne(fname,index)





fileID = fopen(fname);
[outfilpath,outfilname] = fileparts(fname);
fnameout = [outfilpath filesep outfilname '_' num2str(index) '.mat'];
v = fscanf(fileID , '%d');
x = v(index);
v = x + 1;

save(fnameout,'v')
% load fname

% grab the Nth value where N is the index

% add one

%save to disk

end
