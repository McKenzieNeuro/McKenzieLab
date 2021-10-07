function plotMeanSEM(x,y,col,varargin)

if ~isempty(varargin)
    fcn = varargin{1};
else
    fcn = @nanmean;
end
y_m = fcn(y,1);

z = SEM(y,1);
z(isnan(z)) = 0;
y_m(isnan(y_m)) = 0;

plot(x,y_m,'color',col,'linewidth',2)
hold on
plotShadedError(x,y_m,z,col,.5)
hold on



end