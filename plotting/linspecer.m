%% lineStyles=linspecer(N)
% This function creates a cell array of N, [R B G] color matricies
% These can be used to plot lots of lines with distinguishable and nice
% looking colors
%
% lineStyles = linspecer(N);  makes n colors for you to use like:
%                                         plot(x,y,'color',linestyles{ii});
% colormap(linspecer); set your colormap to be a really great one with
%                      easily distinguishable colors and a pleasing aesthetic
% 
% lineStyles = linspecer(N,'qualitative'); forces the colors to all be distinguishable
% lineStyles = linspecer(N,'sequential'); forces the colors to vary along a spectrum 
% lineStyles = linspecer(N,colormap); picks the colors according to your favorite colormap, like 'jet'
% 
%% credits and where the function came from
% The colors are largely taken from:
% http://colorbrewer2.org and Cynthia Brewer, Mark Harrower and The Pennsylvania State University
% 
% She studied this from a phsychometric perspective and crafted the colors
% beautifully.
% 
% I made choices from the many there to decide the nicest once for plotting
% lines in Matlab. I also made a small change to one of the colors I
% thought was a bit too bright. In addition some interpolation is going on
% for the sequential line styles.
% 
%% Example demonstrating the colors.
% C = linspecer(5);
% x = linspace(0,pi*3,1000);
% lgnd = cell(N,1); hold('off');
% for ii=1:N
%     plot(x,sin(x+2*ii*pi/N),'linewidth',5,'color',C{ii}/1);
%     hold('on');
%     lgnd{ii}=num2str(ii);
% end
% legend(lgnd); ylim([-1.1 1.1]);
%% example demonstrating the benefits of the colormap:
% A = rand(15);
% figure; imagesc(A); % default
% figure; imagesc(A); colormap(linspecer); % colorbrewer colormap
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Jonathan Lansey, March 2009-2013  Lansey at gmail.com               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lineStyles=linspecer(N,varargin)

if ~exist('N') % better return a colormap
    temp = linspecer(255);
    temp = [temp{:}];
    lineStyles = reshape(temp,3,255)';
    return;
end

if N<=0 % its empty, nothing else to do here
    lineStyles={};
    return;
end

% interperet varagin
qualFlag = 0;

if ~isempty(varargin)>0 % you set a parameter?
    switch lower(varargin{1})
        case {'qualitative','qua'}
            if N>12 % go home, you just can't get this.
                warning('qualitiative is not possible for greater than 12 items, please reconsider');
            else
                if N>9
                    warning(['Default may be nicer for ' num2str(N) ' for clearer colors use: whitebg(''black''); ']);
                end
            end
            qualFlag = 1;
        case {'sequential','seq'}
            lineStyles = cmap2linspecer(colorm(N));
            return;
        case 'colormap'
                A = colormap;
                if size(A,1)<N
                    warning('your colormap needs more values, there will be duplicate colors');
                end
                values=round(linspace(1,size(A,1),N)); % pick evenly space points along the colormap;
                C = A(values,:); % pick only the colors we want
                lineStyles = cmap2linspecer(C);
                return;
        otherwise
            if exist(varargin{1}) % this is an existing colormap specification. We'll just use that.
                C = eval([varargin{1} '(' num2str(N) ')']);
                lineStyles = cmap2linspecer(C);
                return;
            else
                warning(['parameter ''' varargin{1} ''' not recognized']);
            end
    end
end      
      
% predefine some colormaps
set3 = colorBrew2mat({[141, 211, 199];[ 255, 237, 111];[ 190, 186, 218];[ 251, 128, 114];[ 128, 177, 211];[ 253, 180, 98];[ 179, 222, 105];[ 188, 128, 189];[ 217, 217, 217];[ 252, 205, 229];[ 204, 235, 197];[ 255, 255, 179]}');
set1JL = brighten(colorBrew2mat({[228, 26, 28];[ 55, 126, 184];[ 77, 175, 74];[ 255, 127, 0];[ 255, 237, 111]*.95;[ 166, 86, 40];[ 247, 129, 191];[ 153, 153, 153];[ 152, 78, 163]}'));
set1 = brighten(colorBrew2mat({[ 55, 126, 184]*.95;[228, 26, 28];[ 77, 175, 74];[ 152, 78, 163];[ 255, 127, 0]}),.8);

switch N
    case 1
        lineStyles = { [  55, 126, 184]/255};
    case {2, 3, 4, 5 }
        lineStyles = set1(1:N);
    case {6 , 7, 8, 9}
        lineStyles = set1JL(1:N)';
    case {10, 11, 12}
        if qualFlag % force qualitative graphs
            lineStyles = set3(1:N)';
        else % 10 is a good number to start with the sequential ones.
            lineStyles = cmap2linspecer(colorm(N));
        end
        
otherwise % any old case where I need a quick job done.
    
    lineStyles = cmap2linspecer(colorm(N));
    
end

end
% extra functions
function varIn = colorBrew2mat(varIn)
for ii=1:length(varIn) % just divide by 255
    varIn{ii}=varIn{ii}/255;
end        
end

function varIn = brighten(varIn,varargin) % increase the brightness

if isempty(varargin), frac = .9; else, frac = varargin{1}; end

for ii=1:length(varIn)
    varIn{ii}=varIn{ii}*frac+(1-frac);
end        
end

function vOut = cmap2linspecer(vIn) % changes the format from a double array to a cell array with the right format
vOut = cell(size(vIn,1),1);
for ii=1:size(vIn,1)
    vOut{ii} = vIn(ii,:);
end
end
%%
% colorm returns a colormap which is really good for creating informative
% heatmap style figures.
% No particular color stands out and it doesn't do too badly for colorblind people either.
% It works by interpolating the data from the
% 'spectral' setting on http://colorbrewer2.org/ set to 11 colors
% 
% It is modified a little to make the brightest yellow a little less bright.
% 
% Jonathan Lansey 2013
% % % % % % % % % % % % % % % % 
% % Example:
% x=rand(15);figure;imagesc(x);figure;imagesc(x);colormap(colorm);
% % % 
function cmap = colorm(varargin)
n = 100;
if ~isempty(varargin)
    n = varargin{1};
end

if n==1
    cmap =  [0.2005    0.5593    0.7380];
    return;
end
if n==2
     cmap =  [0.2005    0.5593    0.7380;
              0.9684    0.4799    0.2723];
          return;
end

frac=.95; % Slight modification from colorbrewer here to make the yellows in the center just a bit darker
cmapp = [158, 1, 66; 213, 62, 79; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255*frac, 255*frac, 191*frac; 230, 245, 152; 171, 221, 164; 102, 194, 165; 50, 136, 189; 94, 79, 162];
x = linspace(1,n,size(cmapp,1));
xi = 1:n;
for ii=1:3
    cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
end
cmap = flipud(cmap/255);
end