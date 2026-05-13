
function contourf_alpha(pts, col, alpha_val, n_levels)
    if size(pts,1) < 10, return; end
    x=pts(:,1); y=pts(:,2);
    xi=linspace(min(x)-range(x)*0.2, max(x)+range(x)*0.2, 60);
    yi=linspace(min(y)-range(y)*0.2, max(y)+range(y)*0.2, 60);
    [Xi,Yi]=meshgrid(xi,yi);
    try
        pd = fitdist([x y],'mvn'); % not standard
    catch
    end
    % Manual KDE
    bw = std([x;y])*size(pts,1)^(-1/6);
    Zi = zeros(size(Xi));
    for k=1:size(pts,1)
        Zi = Zi + exp(-((Xi-x(k)).^2+(Yi-y(k)).^2)/(2*bw^2));
    end
    Zi = Zi/sum(Zi(:));
    levels = linspace(prctile(Zi(:),60),max(Zi(:)),n_levels+1);
    [C,h]=contourf(Xi,Yi,Zi,levels,'EdgeColor','none');
    % Apply color and alpha
    kids=get(h,'Children');
    if ~iscell(kids), kids={kids}; end
    for k=1:length(kids)
        try
            frac=k/length(kids);
            set(kids{k},'FaceColor',col,'FaceAlpha',alpha_val*frac,...
                'EdgeColor','none')
        catch, end
    end
end
