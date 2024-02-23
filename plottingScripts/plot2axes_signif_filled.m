function plot2axes_signif_filled(xVals, hVals, pVals)
% plots grey transparent bands of significant values
% xVals = x-axis (e.g. time)
% hVals = 0 (n.s.) or 1 (was significant)
% pVals = p-values

% (c) Jiri, Oct22

hVals = hVals(:);   % column vector N x 1
pVals = pVals(:);   % column vector N x 1

% xVals = trialsData.xVals;
% hVals = trialsData.hVals(i_x,ch);
% pVals = trialsData.pVals(i_x,ch);
% assert(size(hVals,1) == size(trialsData.yVals,1));   % ~ lags
% assert(size(hVals,2) == size(trialsData.yVals,2));   % ~ chnl

yLims_ch = get(gca, 'ylim');
i_signif = find(hVals == 1);
if ~isempty(i_signif)
    for s = i_signif'
        x_delta = mean(diff(xVals));
        xm = xVals(s) - x_delta/2;
        xp = xVals(s) + x_delta/2;
        ym = yLims_ch(1);
        yp = yLims_ch(2);
        C = [0.4, 0.4, 0.4];
%         if pVals(s) < 0.001, C = [0.2 0.2 0.2]; end
        h = fill([xm,xm,xp,xp],[ym,yp,yp,ym], C, 'edgecolor','none');
        alpha(h, 0.3);        
    end
end

% control (uncomment)
% for t = 1:length(i_signif)
%     plot(trialsData.xVals(i_signif(t)), yLims_ch(1), 'k*');
% end
