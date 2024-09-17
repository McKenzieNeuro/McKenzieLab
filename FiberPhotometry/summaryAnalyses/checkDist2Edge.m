
for i = 1:length(dirs)
    
    cd(dirs{i})
    
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        
        con = unique(sessiondata.behavior.position.context);
        [~,b] = histc(sessiondata.behavior.totalEdgDist(:,1),-10:10);
        col = linspecer(21,'jet');
        for j = 1:length(con)
            
        kp = contains(sessiondata.behavior.position.context,con{j});
        kp = kp(1:length(b));
        figure
        for i  =0:20
            plot(sessiondata.behavior.position.left_ear_cor(kp&b==i,1),sessiondata.behavior.position.left_ear_cor(kp&b==i,2),'.','color',col{i+1})
            hold on
        end
        waitforbuttonpress
        close all
        
    end
    
    end
end