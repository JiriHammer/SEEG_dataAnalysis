function plot2axes_paraTimes(paraTimes, clrs, fillBand)
% plots into current axes vertical bands of paradigm times (RT, ...)
% paraTimes = 3D: clz x responseType x "channels" (or repetitions - channels in AA)
% clrs = colors
% if fillBand = 1, plots band = mean +/- SEM  (over repetitions)
% if fillBand = 0, plots vertical line = mean (over repetitions)

% (c) Jiri, Oct22

if nargin < 3
    fillBand = true;
end

% do not plot if out of xlims
xLims = get(gca, 'xlim');
pt_all = mean(paraTimes,3);
pt_all = pt_all(:);
if any(pt_all < xLims(1))
    return;
elseif any(pt_all > xLims(2))
    return;
end

assert(size(paraTimes,1) <= size(clrs,1));
yLims = get(gca, 'ylim');

% plot paraTimes (mean +/- SEM over "channels"/repetitions)
for clz = 1:size(paraTimes,1)
    for tt = 1:size(paraTimes,2)
        t_avg = mean(paraTimes(clz,tt,:),3);
        t_sem = sem(paraTimes(clz,tt,:),3);
        if ~isnan(t_avg)      
            if fillBand
                xm = t_avg - t_sem;
                xp = t_avg + t_sem;
                ym = yLims(1);
                yp = yLims(2);
                C = clrs(clz,:);
                h = fill([xm,xm,xp,xp],[ym,yp,yp,ym], C, 'edgecolor','none');
                alpha(h, 0.3); 
            end
            plot([t_avg, t_avg], yLims, 'Color',clrs(clz,:));
        end
    end
end
