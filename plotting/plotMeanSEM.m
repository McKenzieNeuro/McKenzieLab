function plotMeanSEM(x,y,col,varargin)
% plots the Mean  + shaded errors bars around the mean
% x = Nx1 vector describing the columns of y
% y = M x N matrix where M is the number of observations
% col = 3 digit color spec
%
% OPTIONAL
% yAxisFunction =  function to compute central tendency over rows of y (default = nanmean)
% ErrorBarFunction =  function to compute error bars over rows of y (default = SEM)

p = inputParser;


addParameter(p,'yAxisFunction',@nanmean) 
addParameter(p,'ErrorBarFunction',@SEM) 


parse(p,varargin{:})


ErrorBarFunction = p.Results.ErrorBarFunction;
yAxisFunction = p.Results.yAxisFunction;




y_m = yAxisFunction(y,1);
z = ErrorBarFunction(y,1);

%z = prctile(y,5);
z(isnan(z)) = 0;
y_m(isnan(y_m)) = 0;

plot(x,y_m,'color',col,'linewidth',2)
hold on
plotShadedError(x,y_m,z,col,.5)
hold on



end