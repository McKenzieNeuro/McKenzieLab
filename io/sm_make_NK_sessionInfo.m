function sm_make_NK_sessionInfo(fname)


outfil = [fname(1:end-3) 'mat'];
[sFile, ChannelMat] = in_fopen_nk(fname);

sessionInfo.sFile = sFile;
sessionInfo.ChannelMat = ChannelMat;


save(outfil,'sessionInfo')
end