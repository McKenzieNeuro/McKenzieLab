function X= upSample (input, bins)
bins=round(bins);
input=input(:);
X=nan(length(input)*bins,1);


for i = 0:bins-1
X(1+i:bins:end)=input;
end


end