function info = viz(K, varargin)
%VIZ Visualize a ktensor.
%
%   VIZ(K) visualizes the components of an D-way ktensor with R components
%   in an R x D arrangment of plots. Each column of plots represents the
%   columns of the associated factor matrix, and each row represents a
%   single rank-one component.
%
%   INFO = VIZ(K, parameter, value...) takes optional parameters and
%   returns additional information from the plot, including handles to all
%   the axes.
%
%   Primary optional parameters:
%
%   'PlotCommands'- Plot command, 'plot','bar','scatter', or a function
%                   handle of the form |@(x,y) plot(x,y);|.  
%                   Can also be a cell array, one entry per mode.
%   'ModeTitles' - Cell array of mode titles. Use 'none' to disable.
%   'Figure' - Figure number. Default: [] (new figure).
%   'RelModeWidth' - Relative vertical width per mode. Default: ones(D,1).
%   'FactorTitles' - Choices are 'Weight' (relative), 'True Weight',
%                    'Number' or 'None'. Default: 'Weight'.
%   'Normalize' - Automatically normalizes the factor matrix columns and
%                 sorts the components by weight. Options: -1 (none), 1
%                 (1-norm), 2 (2-norm). Can also be a function handle to a
%                 function that accepts a ktensor as input. Default: 2. 
%   'Title' - Title for the overall plot. If used together with mode
%             titles, then adjust 'TopSpace' accordingly.
%
%   Detailed optional parameters:
%
%   -- Spacing (all proportions in [0,1]) --
%   'TopSpace' - Space at the top, for titles. Default: 0.05.
%   'BottomSpace' - Space at bottom, for xticklabels. Default: 0.05.
%   'LeftSpace' - Space at left, for component labels. Default: 0.025.
%   'RightSpace' - Space at right. Default: 0.025.
%   'VertSpace' - Vertical space inbetween factor axes. Default: 0.01.
%   'HorzSpace' - Horizontal space inbetween factor axes. Default: 0.01.
%   'XVals' - X-values, specify one per mode:
%             o [] - Default 1:size(K,n)
%             o xvals - Array of length size(K,n).
%             Default: repmat({[]},[ndims(K) 1]).
%   'XLims' - Choose one per mode:
%             o 'same' - Same x-limits on all axes in the same mode
%             o [xmin xmax] - User-specified x-limits for all axes
%             o [] - No modification to what is done by the plot routine
%             Default: repmat({'same'},[ndims(K) 1]).
%   'YLims' - Choose one per mode:
%             o 'same' - Same y-limits on all axes (including zero)
%             o [ymin ymax] - User-specified y-limits for all axes
%             o [] - No modification to what is done by the plot routine
%             o 'addzero' - Different y-limits but adjust so zero shown
%             Default: repmat({'same'},[ndims(K) 1]).
%   'ShowZero' - Show dashed line at zero? Default: true(ndims(K),1).
%   'XTicks' - Boolean array for showing xticks. Default: true(ndims(K),1).
%   'YTicks' - Boolean array for showing yticks. Default: false(ndims(K),1).
%              (If YTicks is true, then need to increase 'HorzSpace'.)
%   'BaseFontSize' - Smallest font size. Default: 14.
%
%   Return values:
%   'height' - Height of each plot (as a proportion in [0,1]).
%   'width' - Width of each plot (as a proportion in [0,1]).
%   'ModeTitles' - Handles for the D mode titles.
%   'GlobalAxis' - Handle for main axes in which all others are embedded.
%   'FactorAxes' - D x R array of handles for the subfigure axes.
%   'ModeTitleHandles' - D-array of the handles to the mode titles (on the top).
%   'CompTitleHandles' - R-array of handles to the factor titles (on the left).
%   'PlotHandles'- D x R array of handles to the figures for each factor.
%
%   Examples:
%   K = ktensor([3; 2], rand(40,2), rand(50,2), rand(30,2));
%   viz(K,'Figure',1,'Hspace',0.05,'Vspacebottom',0.075);
%
%   Thanks to Alex Williams for the prototype for this functionality.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% TGK: Need to add options around line at zero, marks to use, font sizes, etc.



%%
nd = ndims(K); % Order
nc = ncomponents(K); % Rank

%% Parse optional inputs

params = inputParser;
% Normalize
params.addParameter('Normalize',2);
% Figure
params.addParameter('Figure', []);
% Spacing
params.addParameter('RelModeWidth', ones(nd,1)); % Horizontal space for each mode
params.addParameter('HorzSpace',0.01); % Horizontal space between axes
params.addParameter('RightSpace',0.025); % Horizontal space on left
params.addParameter('LeftSpace',0.075); % Horizontal space on right
params.addParameter('VertSpace',0.01); % Vertical space between axes
params.addParameter('TopSpace',0.05); % Vertical space at top
params.addParameter('BottomSpace',0.05); % Vertical space at bottom
% Titles
params.addParameter('ModeTitles', 'default');
params.addParameter('FactorTitles', 'weight'); % Default is 'none'. Options are 'weight' or 'number'
params.addParameter('Title','none'); % Overall title
% Plots
params.addParameter('PlotCommands', ...
    repmat({@(x,y) plot(x,y,'LineWidth',1,'Color','b');}, [nd 1]));
params.addParameter('XVals',repmat({[]},[nd 1]));
params.addParameter('XLims', repmat({'same'},[nd 1]));
params.addParameter('YLims', repmat({'same'},[nd 1]));
params.addParameter('ShowZero',true(nd,1));
params.addParameter('XTicks',true(nd,1));
params.addParameter('YTicks',false(nd,1));
params.addParameter('BaseFontSize',14);


params.parse(varargin{:});
res = params.Results;

%% Handle scalar 'XTicks', 'YTicks'
if isscalar(res.XTicks)
   res.XTicks = repmat(res.XTicks,[nd 1]);
end
if isscalar(res.YTicks)
   res.YTicks = repmat(res.YTicks,[nd 1]);
end

%% Clean up tensor - applying normalization
if isa(res.Normalize,'function_handle')
    K = res.Normalize(K);
elseif res.Normalize > 0
    fprintf(['ktensor/viz: Normalizing factors and sorting components ' ...
        'according to the %d-norm.\n'], res.Normalize);
    K = normalize(K,'sort',res.Normalize);
end

%% Create new figure or reset old figure
if isempty(res.Figure)
    figure;
else
    figure(res.Figure);
    clf;
end

%% Create number-of-modes (nd) x number-of-components (nc) axes

% Calculate the amount of vertical space available for the plots themselves
% by subtracting off the top and bottom space as well as the inbetween
% space.
Vplotspace = 1 - res.TopSpace - res.BottomSpace - (nc - 1) * res.VertSpace;
height = Vplotspace / nc;

% Do likewise for the horizontal space.
Hplotspace = 1 - res.LeftSpace - res.RightSpace - (nd - 1) * res.HorzSpace;
width = (res.RelModeWidth ./ sum(res.RelModeWidth)) .* Hplotspace;

% Create the global axis
GlobalAxis = axes('Position',[0 0 1 1]); % Global Axes
axis off;

% Create the nd x nc factor axes array
FactorAxes = gobjects(nd,nc); % Factor Axes
for k = 1 : nd
    for j = 1 : nc
        xpos = res.LeftSpace + (k-1) * res.HorzSpace + sum(width(1:k-1));
        ypos = 1 - res.TopSpace - height - (j-1) * (height + res.VertSpace);
        FactorAxes(k,j) = axes('Position',[xpos ypos width(k) height]);
        FactorAxes(k,j).FontSize = res.BaseFontSize;
    end
end


%% Plot each factor

% Set up the plot command for each mode
PlotCommands = res.PlotCommands;
if ~iscell(PlotCommands)
    PlotCommands = repmat({PlotCommands}, [nd 1]);
end
for k = 1:nd
    if isempty(PlotCommands{k}) || strcmpi(PlotCommands{k},'plot') ...
            || strcmpi(PlotCommands{k},'line')
        PlotCommands{k} = @(x,y) plot(x,y,'Linewidth',1,'Color','b');
    elseif strcmpi(PlotCommands{k},'bar')
        PlotCommands{k} = @(x,y) bar(x,y,'EdgeColor','b','FaceColor','b');
    elseif strcmpi(PlotCommands{k},'scatter')
        PlotCommands{k} = @(x,y) scatter(x,y,10,'b','filled');
    end
end

% Set up cell array to store handles for all the plots
h = cell(nd,nc);

% Loop through each mode 
for k = 1 : nd

    % Extract kth factor matrix from K
    U = K.u{k};

    % Create x-data for this mode
    if isempty(res.XVals{k})
        xx = (1:size(K,k))';
    else
        xx = res.XVals{k};
        % Check that xx is in strictly increading order
        if any(diff(xx) <= 0)
            error('XVals for mode %d are not in strictly increasing order', k)
        end
    end

    % Adjust default xlims for this mode
    if (xx(end) - xx(1)) >= 2
        xdelta = [-1 1];
    else
        tmp = 0.1 *(xx(end) - xx(1));
        xdelta = [-tmp tmp];
        clear tmp
    end
    xl = [xx(1) xx(end)] + xdelta;

    % Create default ylims for this mode
    if isequal(res.YLims{k},'addzero')
        yl = [min( -0.01, min(U(:)) ), max( 0.01, max(U(:)) )];
    else
        yl = [min( 0, min(U(:)) ), max(U(:)) ];
    end
    % Loop through each factor in mode
    for j = 1 : nc

        % Extract y data
        yy = U(:,j);

        % Set up plot
        hold(FactorAxes(k,j), 'off');
        axes(FactorAxes(k,j));

        % Call the plot command
        hh = PlotCommands{k}(xx, yy);

        % Set x-axes
        if isequal(res.XLims{k}, 'same')
            xlim(FactorAxes(k,j),xl);
        elseif isnumeric(res.XLims{k}) && isequal(size(res.XLims{k}),[1 2])
            xlim(FactorAxes(k,j),res.XLims{k});
        else
            fprintf(['ktensor/viz: Doing nothing to FactorAxes xlim' ...
                ' for mode %d and component %d\n'],k,j);
        end

        % Set y-axes
        if isequal(res.YLims{k}, 'same')
            ylim(FactorAxes(k,j),yl);
        elseif isequal(res.YLims{k},'addzero')
            % Create y-axes that include zero
            tmpyl = [ min(-0.01, min(U(:,j))), max( 0.01, max(U(:,j))) ];
            ylim(FactorAxes(k,j),tmpyl);
        elseif isnumeric(res.YLims{k}) && isequal(size(res.YLims{k}),[1 2])
            ylim(FactorAxes(k,j),res.YLims{k});
        else
            fprintf(['ktensor/viz: Doing nothing to FactorAxes ylim' ...
                ' for mode %d and component %d\n'],k,j);
        end

        % Turn off y-label
        set(FactorAxes(k,j), 'Ylabel', []);

         % Turn off x-ticks?
        if ~res.XTicks(k)
            set(FactorAxes(k,j),'Xtick',[]);
        end

        % Turn off y-ticks?
        if ~res.YTicks(k)
            set(FactorAxes(k,j),'Ytick',[]);
        end

        % Draw a box around the axes
        set(FactorAxes(k,j),'Box','on')

        % Turn off X ticklabels if not the bottom plot
        if j < nc
            set(FactorAxes(k,j),'XtickLabel',{});
        end

        % Draw dashed line at zero?
        if res.ShowZero(k)
            xl = xlim(FactorAxes(k,j));
            hold(FactorAxes(k,j), 'on');
            plot(FactorAxes(k,j), xl, [0 0], 'k:', 'Linewidth', 1.5);
        end

        % Save handle for main plot
        h{k,j} = hh;

        % Make the fonts on the xtick labels big
        set(FactorAxes(k,j),'FontSize',res.BaseFontSize)
    end
end

%% Title for each mode, along the top
if ( isa(res.ModeTitles,'char') && strcmpi(res.ModeTitles,'none') )
    % No titles
    ModeTitleHandles = repmat({[]},nd,1);
else
    % Create handles for the titles
    ModeTitleHandles = gobjects(nd,1);

    % Set ModeTitles if needed
    if ( isa(res.ModeTitles,'char') && strcmpi(res.ModeTitles,'default') )
        ModeTitles = cell(nd,1);
        for i = 1:nd
            ModeTitles{i} = sprintf('Mode %d',i);
        end
    else
        ModeTitles = res.ModeTitles;
    end

    % Add the mode titles
    axes(GlobalAxis);
    for k = 1:nd
        xpos = res.LeftSpace ...
            + (k-1) * res.HorzSpace + sum(width(1:k-1)) ...
            + 0.5 * width(k);
        ypos = 1 - res.TopSpace;
        ModeTitleHandles(k) = text(xpos,ypos,ModeTitles{k},...
            'VerticalAlignment','Bottom','HorizontalAlignment','Center');
        set(ModeTitleHandles(k),'FontSize',res.BaseFontSize+1)
        set(ModeTitleHandles(k),'FontWeight','bold')
    end
end

%% Print component titles along the left side
CompTitleHandles = gobjects(nc,1);
if ~strcmpi(res.FactorTitles,'none')
    axes(GlobalAxis);
    rellambda = abs (K.lambda / K.lambda(1));
    for j = 1:nc
        xpos = 0.1 * res.LeftSpace;
        ypos = 1 - res.TopSpace - 0.5 * height ...
            - (j-1) * (height + res.VertSpace);
        if strcmpi(res.FactorTitles,'weight')
            txt = sprintf('%3.2f', rellambda(j));
        elseif strcmpi(res.FactorTitles,'true weight')
            txt = sprintf('%3.2f', K.lambda(j));
        else
            txt = sprintf('%d', j);
        end
        CompTitleHandles(j) = text(xpos,ypos,txt,...
            'VerticalAlignment','Middle','HorizontalAlignment','Left');
        set(CompTitleHandles(j),'FontSize',res.BaseFontSize)
    end
end

%% Global Title
TitleHandle = gobjects(1,1);
if ~strcmpi(res.Title,'none')
    axes(GlobalAxis);
    TitleHandle = text(0.5,1,res.Title,...
        'VerticalAlignment','Top','HorizontalAlignment','Center');
    set(TitleHandle,'FontSize',res.BaseFontSize+2)
    set(TitleHandle,'FontWeight','bold')
end


%% Save stuff to return
info.height = height;
info.width = width;
info.GlobalAxis = GlobalAxis;
info.FactorAxes = FactorAxes;
info.ModeTitleHandles = ModeTitleHandles;
info.CompTitleHandles = CompTitleHandles;
info.PlotHandles = h;
info.TitleHandle = TitleHandle;
