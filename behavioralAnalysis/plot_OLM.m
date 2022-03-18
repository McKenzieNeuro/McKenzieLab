% get all mat files

masterDir = 'R:\DonaldsonT\OLM-210922-100555';

fils = getAllExtFiles(masterDir,'mat',1);
fils = fils([3 6:end]);
for i = 1:length(fils)
    
   d = load(fils{i}); 
   % get object 1 on/off
   
   ix1 = cellfun(@any,regexpi(d.data(:,1),'Obj1on'));
   ix2 = cellfun(@any,regexpi(d.data(:,1),'Obj1off'));
   ix3 = cellfun(@any,regexpi(d.data(:,1),'Obj2on'));
   ix4 = cellfun(@any,regexpi(d.data(:,1),'Obj2off'));
   
   Obj1on = cell2mat(d.data(ix1,2));
   Obj1off = cell2mat(d.data(ix2,2));
   Obj2on = cell2mat(d.data(ix3,2));
   Obj2off = cell2mat(d.data(ix4,2));
   
   Obj1 = [Obj1off  - Obj1on];
   Obj2 = [Obj2off -  Obj2on];
   
   [max(Obj1) min(Obj1)]
   
   [max(Obj2) min(Obj2)]
   
   
   d_prime(i) = (sum(Obj2) - sum(Obj1)) / (sum(Obj2) + sum(Obj1));
   
   i
end