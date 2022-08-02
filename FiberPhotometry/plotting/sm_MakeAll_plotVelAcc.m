fils = getAllExtFiles('R:\DANEHippocampalResponse','csv',1);
filenameps = 'R:\DANEHippocampalResponse\Figures\allDANE_vel_acc.ps';
filenamepdf = 'R:\DANEHippocampalResponse\Figures\allDANE_vel_acc.pdf';

[dirs,b] = fileparts(fils);
dirs = unique(dirs);

for i = 1:length(dirs)
    try
        h = sm_plotVelAcc(dirs{i});
        print(h, '-dpsc2',filenameps ,'-append')
        close all
    catch
        disp('bad dir')

    end
    i
end

%%

ps2pdf('psfile',filenameps,'pdffile',filenamepdf)

