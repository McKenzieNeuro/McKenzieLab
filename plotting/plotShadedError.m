function plotShadedError(x,y,z,col,alpha)


x = x(:);
y = y(:);
z = z(:);
y(isnan(y)) = 0;
z(isnan(z)) = 0;
xx = [x;flipud(x)];
zz = [y-z;flipud(y)+flipud(z)];



patch(xx,zz,col,'facealpha',alpha,'edgecolor','none')


end