function plotTrials_yy(params, trialsData)
% double y-axis plots, plotyy mode (2 y-axis)
% plots values in 'trialsData' struct with fields:
%    yVals_1: [sample x ch]
%    yErrs_1: [sample x ch x clz] -> optional (plotband)
%    xVals_1: x-Axis values (vector), typically time in [s]
% and
%    yVals_2: [sample x ch]
%    yErrs_2: [sample x ch x clz] -> optional (plotband)
%    xVals_2: x-Axis values (vector), typically time in [s]
%    info: struct with supplementary information

% (c) Jiri, May17

%% basic check
nCh = size(trialsData.yVals_1,2);
assert(size(trialsData.yVals_1,2) == size(trialsData.yVals_2,2));
assert(size(trialsData.yVals_1,1) == size(trialsData.xVals_1,1));
assert(size(trialsData.yVals_2,1) == size(trialsData.xVals_2,1));

%% load 'H'
%load(params.storage.cacheFile, 'H');

%% y-lim
% if isfield(trialsData.info, 'yLims')
%     yLims = trialsData.info.yLims;      % if empty, no common y-limits!
% else
%     if isfield(trialsData, 'yErrs')
%         yLims_lo = getYLims(trialsData.yVals - trialsData.yErrs);
%         yLims_hi = getYLims(trialsData.yVals + trialsData.yErrs);
%         yLims = [yLims_lo(1), yLims_hi(2)];
%     else
%         yLims = getYLims(trialsData.yVals);
%     end
% end

%% font size
if isfield(trialsData.info, 'fontsize')
    fontsize = trialsData.info.fontsize;
else
    fontsize = 12;
end

%% outliers
% if isfield(trialsData.info, 'outliers')
%     outliers = trialsData.info.outliers;
% else
%     outliers = getOutliers(trialsData.yVals);
% end
% assert(length(outliers) == nCh);
%outliers = zeros(nCh,1);    % hack!

%% figure
%f = figure('visible','on', 'Position',params.plot_general.scrnsz_pixels);
%f = figure('visible','on', 'Position',[1 -299 1920 932]);
f = figure('visible','on', 'Position',[1 -479 2880 1463]);

[nRows, nCols] = getSubplotLayout(nCh);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.005, 0.005];

%% text on upper part of the figure
if isfield(trialsData.info, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = trialsData.info.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fonts', 16, 'fontw', 'bold');
end

%% colors
% if isfield(trialsData.info, 'colors')
%     colors = trialsData.info.colors;
% else
%     colors = colorPalette(size(trialsData.yVals,3));
% end

%% plot channel-wise
for ch = 1:nCh

    % if outlier, set background color to pink
%     if any(outliers(ch,:))
%         bkgClr = [1 0.8 0.8];
%         txtClr = 'r';
%     else
%         bkgClr = [1 1 1];
%         txtClr = 'k';
%     end
    txtClr = 'k';
    bkgClr = [1 1 1];
    isOutlier = false;
    
    % new subplot   
    %subplot(nRows, nCols, ch);
    subtightplot(nRows, nCols, ch, gap, marg_h, marg_w);
    set(gca, 'Color',bkgClr);
    set(gca, 'FontSize',fontsize);
    hold on;

    % >>>>>>>>>>>>>>> plot <<<<<<<<<<<<<<<<<<<<
    [hAx,hLine1,hLine2] = plotyy(trialsData.xVals_1, trialsData.yVals_1(:,ch),...
                                 trialsData.xVals_2, trialsData.yVals_2(:,ch));
    
   
    % set x-limits
    set(hAx(1), 'xlim',[trialsData.xVals_1(1), trialsData.xVals_1(end)]);
    set(hAx(2), 'xlim',[trialsData.xVals_2(1), trialsData.xVals_2(end)]);
    set(hAx(2), 'XTick',[]);
    xLims = [trialsData.xVals_2(1), trialsData.xVals_2(end)];
    
    % set y-limits
    set(hAx(1), 'ylim',trialsData.info.yLims_1);
    yLims_ch = (get(hAx(1),'ylim'));
    if max(trialsData.yVals_2(:,ch)) <= trialsData.info.yLims_2(end)
        set(hAx(2), 'ylim',trialsData.info.yLims_2);
    else
        isOutlier = true;
        txtClr = 'r';
    end
    
    % click on plot to see it bigger
    set(hAx(1), 'ButtonDownFcn','call_copy');
    set(hAx(2), 'ButtonDownFcn','call_copy');
    
    % channel name
    if isfield(trialsData.info, 'chNames')
        chName = trialsData.info.chNames{ch};
    else
        chName = ['channel: ' num2str(ch)];
    end
    %chName = strrep(chName, '_','\_');
    %title(chName, 'FontSize',9, 'FontWeight','bold');
    pos_text = [xLims(1)+0.5*diff(xLims), yLims_ch(1)+0.9*diff(yLims_ch)];
    text(pos_text(1),pos_text(2), chName, ...
        'FontSize',fontsize, 'FontWeight','bold', 'HorizontalAlignment','center', 'Color',txtClr, 'Interpreter','none');     % 
    
    % xlabel (bottom line)
    if ch > (nRows-1)*nCols
        if isfield(trialsData.info, 'xlabel')
            xlabel(trialsData.info.xlabel, 'FontSize',fontsize);
        end
    else
        set(gca, 'XTick',[]);
    end

    % ylabel (left column)
    if mod(ch,nCols) == 1 || nCh == 1
        if isfield(trialsData.info, 'ylabel')
            ylabel(trialsData.info.ylabel_1, 'FontSize',fontsize);
        end
    elseif mod(ch,nCols) == 0 || nCh == 1
        if isfield(trialsData.info, 'ylabel')
            ylabel(trialsData.info.ylabel_2, 'FontSize',fontsize);
        end
    else
        set(hAx(1), 'YTick',[]);
        if ~isOutlier
            set(hAx(2), 'YTick',[]);
        end
    end
    
    % vertical & horizontal bars (if specified)
    if isfield(trialsData.info, 'verLines')
        for t = trialsData.info.verLines
            plot(hAx(1), [t, t], yLims_ch, '--b');
        end
    end
    if isfield(trialsData.info, 'horLines')
        for t = trialsData.info.horLines
            plot([xLims(1), xLims(2)], [t, t], '--k');
        end
    end        
end

%% output directory
if isfield(trialsData.info, 'outDir')
    outDir = trialsData.info.outDir;
else
    outDir = params.storage.outputDir;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

%% figure name
if isfield(trialsData.info, 'figName')
    figname = trialsData.info.figName;
else
    figname = 'notNamed';
end

%% save
set(f, 'PaperPositionMode','auto');
saveas(f, [outDir filesep figname '.fig']);
print(f, '-dpng','-r0', [outDir filesep figname '.png']);
%print(f, '-dpng','-r600', [outDir filesep figname '.png']);
%print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname '.tif']);
close(f);    
display(['Figure: ' figname ' stored in: ' outDir]);

