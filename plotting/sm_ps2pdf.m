function sm_ps2pdf(filenameps,filenamepdf,gscommand)
%print(h, '-dpsc2',filenameps ,'-append');


if isempty(gscommand)
    gscommand = 'C:\Program Files (x86)\gs\gs9.54.0\bin\gswin32.exe';
end
ps2pdf('psfile', filenameps ,'pdffile',filenamepdf,'gscommand',gscommand)


end