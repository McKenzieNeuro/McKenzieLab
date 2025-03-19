function [Out] = symsepchar(StrIn,Sym)

% Takes a string of characters (StrIn) seperated by any special character (Sym)
% and outputs a double array. ex 'C:\Users\Data.txt' -> [C:\,Users\Data.txt] or '11/30/2016' -> [11,30,2016]

symcount = 1;
charcount = 1;
for ci = 1:length(StrIn)
    if strcmp(StrIn(1,ci),Sym) == 1
        symcount = symcount+1;
        charcount = 1;
        continue
    elseif strcmp(StrIn(1,ci),Sym) == 0 
        Out{1,symcount}(1,charcount) = StrIn(1,ci); 
        charcount = charcount+1;
    end
end


