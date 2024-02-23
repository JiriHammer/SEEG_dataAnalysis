function wasOutOfScale = plotChannel2Axes(trialsData, ch, plotInfo)
% plots channel activations (time-series) to specified axes

% (c) Jiri, Apr18

wasOutOfScale = false;      % default (if the channel values were out of yLims scale)

%% settings (defaults)
if ~isfield(plotInfo,'yLims'), yLims = []; else, yLims = plotInfo.yLims; end
if ~isfield(plotInfo,'txtClr'), plotInfo.txtClr = 'k'; end
if ~isfield(plotInfo,'subtractAvgResponse'), plotInfo.subtractAvgResponse = false; end
if ~isfield(plotInfo,'time2plot'), plotInfo.time2plot = [trialsData.xVals(1), trialsData.xVals(end)]; end

%% axes
h = plotInfo.chnlsAxes;

%% colors
if isfield(trialsData.info, 'colors')
    colors = trialsData.info.colors;
else
    colors = colorPalette(size(trialsData.yVals,3));
end

%% font size
if isfield(trialsData.info, 'fontsize')
    fontsize = trialsData.info.fontsize;
else
    fontsize = 12;
end
set(gca, 'FontSize',fontsize);

%% x-range
% if isfield(trialsData.info, 'time2plot')
%     i_x = closestval(trialsData.xVals, trialsData.info.time2plot(1)):closestval(trialsData.xVals, trialsData.info.time2plot(2));
% else
%     i_x = 1:size(trialsData.xVals,1);
% end

%% data -> x, y_avg, y_sem
% trialsData.xVals = trialsData.xVals(i_x);
% trialsData.yVals = trialsData.yVals(i_x,:,:);
% if isfield(trialsData, 'yErrs')
%     trialsData.yErrs = trialsData.yErrs(i_x,:,:);
% end
x = trialsData.xVals;
y_avg = trialsData.yVals(:,ch,:);
if isfield(trialsData, 'yErrs')
    y_sem = trialsData.yErrs(:,ch,:);
end
if plotInfo.subtractAvgResponse
    y_avg = y_avg - mean(y_avg,3);  % subtract AVG (non-specific) response (mean over clz)
end
            
%% >>>>>>>>>>>>>>> plot <<<<<<<<<<<<<<<<<<<<
for clz = 1:size(trialsData.yVals,3)
    clr = colors(clz,:);
%     if isfield(trialsData, 'yErrs')
%         plotband(trialsData.xVals, trialsData.yVals(:,ch,clz), trialsData.yErrs(:,ch,clz), clr, 0.5, 'LineWidth',2);    %  mean - error band plot
%     else
%         plot(trialsData.xVals, trialsData.yVals(:,ch,clz), 'Color',clr, 'LineWidth',0.5);
%     end
    if isfield(trialsData, 'yErrs')
        plotband(x, y_avg(:,1,clz), y_sem(:,1,clz), clr, 0.5, 'LineWidth',2);    %  mean - error band plot
    else
        plot(x, y_avg(:,1,clz), 'Color',clr, 'LineWidth',0.5);
    end    
end

%% plot the AVG (non-specific) response (mean over clz)
if ~plotInfo.subtractAvgResponse
    y_nonSpec = mean(y_avg,3);  % mean over clz
    plot(x, y_nonSpec, ':k', 'LineWidth',0.5);
end

%% set y-limits
axis tight;
xlim(plotInfo.time2plot);
xLims = get(gca, 'xlim');
yLims_ch = yLims;
if ~isempty(yLims_ch)
    tmp = squeeze(trialsData.yVals(:,ch,:));
    if yLims_ch(1) > min(tmp(:)) || yLims_ch(2) < max(tmp(:))   % channel vals out of limits
        plotInfo.txtClr = 'r';
        set(gca, 'ylim', 1.1*get(gca,'ylim')); 
        %outliers(ch,:) = 1;
        wasOutOfScale = true;
    else
        set(gca, 'ylim',yLims_ch); 
    end
end
yLims_ch = double(get(gca, 'ylim'));

%% mark rejected indices (typically for single trial plots)
if isfield(trialsData.info, 'rejected')
    i_reject = find(trialsData.info.rejected(:,ch) == 1);
%     x = trialsData.xVals;
%     ype = yLims_ch(2)*ones(size(x,1),1);
%     yme = yLims_ch(2)*ones(size(x,1),1);
%     yme(i_reject) = yLims_ch(1);
%     h(1,1)= fill([x;x(end:-1:1)],[ype;yme(end:-1:1)], [1.0 0.0 0.0], 'edgecolor','none');
%     alpha(h(1,1), 0.5);        
    for t = 1:length(i_reject)
        plot(trialsData.xVals(i_reject(t)), yLims_ch(1), 'r*');
    end                
end

%% click on plot to see it bigger
set(gca, 'ButtonDownFcn','call_copy');

%% significance between 2 classes
if isfield(trialsData, 'hVals')
    assert(size(trialsData.hVals,2) == size(trialsData.yVals,2));   % ~ chnl 
    assert(size(trialsData.hVals,1) == size(trialsData.yVals,1));   % ~ lags
    hVals = trialsData.hVals(:,ch);
    pVals = trialsData.pVals(:,ch);
    plot2axes_signif_filled(x, hVals, pVals);
    
%     trialsData.hVals = trialsData.hVals(i_x,:);
%     trialsData.pVals = trialsData.pVals(i_x,:);
%     assert(size(trialsData.hVals,1) == size(trialsData.yVals,1));   % ~ lags
%     assert(size(trialsData.hVals,2) == size(trialsData.yVals,2));   % ~ chnl
%     i_signif = find(trialsData.hVals(:,ch) == 1);
%     if ~isempty(i_signif)
%         for s = i_signif'
%             x_delta = mean(diff(trialsData.xVals));
%             xm = trialsData.xVals(s) - x_delta/2;
%             xp = trialsData.xVals(s) + x_delta/2;
%             ym = yLims_ch(1);
%             yp = yLims_ch(2);
%             C = [0.4, 0.4, 0.4];
%             if trialsData.pVals(s,ch) < 0.001, C = [0.2 0.2 0.2]; end
%             h = fill([xm,xm,xp,xp],[ym,yp,yp,ym], C, 'edgecolor','none');
%             alpha(h, 0.3);        
%         end
%     end
%     for t = 1:length(i_signif)
%         plot(trialsData.xVals(i_signif(t)), yLims_ch(1), 'k*');
%     end        
end

%% significance versus baseline
if isfield(trialsData, 'hVals_base')
    assert(size(trialsData.hVals_base,1) == size(trialsData.yVals,1));   % ~ lags
    assert(size(trialsData.hVals_base,2) == size(trialsData.yVals,2));   % ~ chnl
    assert(size(trialsData.hVals_base,3) == size(trialsData.yVals,3));   % ~ clzs
    for clz = 1:size(trialsData.hVals_base,3)
        clr = colors(clz,:);
        i_signif = find(trialsData.hVals_base(:,ch,clz) == 1);
        for t = 1:length(i_signif)
            plot(trialsData.xVals(i_signif(t)), yLims_ch(1)+0.05*clz*diff(yLims_ch), '.', 'MarkerEdgeColor',clr, 'MarkerFaceColor',clr, 'MarkerSize',12);
        end        
    end
end

%% set again xlim
xlim(plotInfo.time2plot);
xLims = get(gca, 'xlim');

%% text (title): channel name
if isfield(trialsData.info, 'chNames')
    chName = trialsData.info.chNames{ch};
else
    chName = ['channel: ' num2str(ch)];
end
if size(chName,2) > 20
    i_semiColon = strfind(chName,':');
    if ~isempty(i_semiColon)
        tmp = chName;
        chName = cell(2,1);
        chName{1} = tmp(1:i_semiColon);
        chName{2} = tmp(i_semiColon+1:end);
    end
end
%chName = strrep(chName, '_','\_');
%title(chName, 'FontSize',9, 'FontWeight','bold');
if ~isempty(chName)
    pos_text = [xLims(1)+0.5*diff(xLims), yLims_ch(1)+0.9*diff(yLims_ch)];
    text(pos_text(1),pos_text(2), chName, ...
        'FontSize',fontsize, 'FontWeight','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr, 'Interpreter','none');     % 
end

%% vertical & horizontal lines (if specified)
if isfield(trialsData.info, 'verLines')
    for t = trialsData.info.verLines
        if t > xLims(1) && t < xLims(2)
            plot([t, t], yLims_ch, '--k');
        end
    end
end
if isfield(trialsData.info, 'horLines')
    for t = trialsData.info.horLines
        plot([xLims(1), xLims(2)], [t, t], '--k');
    end
end        

%% reaction times as vertical lines (if specified)
if isfield(trialsData.info, 'paraTimes')
    for clz = 1:size(trialsData.info.paraTimes,1)
        clr = colors(clz,:);
        for t = trialsData.info.paraTimes(clz,:)
            if ~isnan(t)
                if t > xLims(1) && t < xLims(2)
                    plot([t, t], yLims_ch, 'Color',clr, 'LineWidth',2);
                end
            end
        end
    end
end

