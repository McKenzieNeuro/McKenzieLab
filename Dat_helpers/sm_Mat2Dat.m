function sm_Mat2Dat(mat,fnameOut)

%inputs
%mat = data file, rows = channels
% fnameout = *.dat file with int16 for each channel merged



fidO = fopen(fnameOut,'w');

fwrite(fidO,mat(:),'int16');

fclose(fidO);




end