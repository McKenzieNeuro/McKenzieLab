function Y = nanBuffer(X,siz)


 Y = [nan(size(X,1),siz) X nan(size(X,1),siz)];
 
end