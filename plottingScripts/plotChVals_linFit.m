function slopeFit = plotChVals_linFit(params, dataStruct)
% plots values per channels with a linear fit
% returns slope of linear fits

% (c) Jiri, Sep16

%% basic check
nCh = size(dataStruct.yVals,2);
assert(size(dataStruct.yVals,1) == size(dataStruct.xVals,1));

%% load 'H'
load(params.storage.cacheFile, 'H');

%% y-lim
if isfield(dataStruct.info, 'yLims')
    yLims = dataStruct.info.yLims;
else
    yLims = getYLims(dataStruct.yVals);
end

%% figure
f = figure('visible','on', 'Position',params.plot_general.scrnsz_pixels);

[nRows, nCols] = getSubplotLayout(nCh);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.005, 0.005];

%% text on upper part of the figure
if isfield(dataStruct.info, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = dataStruct.info.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fonts', 14, 'fontw', 'bold');
end

%% colors
if isfield(dataStruct.info, 'colors')
    colors = dataStruct.info.colors;
else
    colors = colorPalette(size(dataStruct.yVals,3));
end
bkgClr = 'w';

%% plot channel-wise
slopeFit = nan(1,nCh);
for ch = 1:nCh

    % new subplot   
    %subplot(nRows, nCols, ch);
    subtightplot(nRows, nCols, ch, gap, marg_h, marg_w);
    set(gca, 'Color',bkgClr);
    hold on;

    % plot
    plot(dataStruct.xVals(:,ch), dataStruct.yVals(:,ch), 'bo');
    
    % axis limits
    set(gca, 'ylim',yLims);
    xLims = get(gca, 'xlim');
    
    % linear fit
    stats = regstats(dataStruct.yVals(:,ch), dataStruct.xVals(:,ch));
    C = stats.beta;
    xFit = [xLims(1):0.01*diff(xLims):xLims(2)];
    yFit = C(1) + C(2)*xFit;
    plot(xFit, yFit, 'Color','r', 'LineWidth',1)
    slopeFit(ch) = C(2);
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn','call_copy');
    
    % significance?
    if isfield(dataStruct, 'hVals')
        i_signif = find(dataStruct.hVals(:,ch) == 1);
        for t = 1:length(i_signif)
            plot(dataStruct.xVals(i_signif(t)), 0.9*max(get(gca,'yLim')), 'r*');
        end
    end

    % channel name
    if isfield(dataStruct.info, 'chNames')
        chName = dataStruct.info.chNames{ch};
    else
        chName = ['channel: ' num2str(ch)];
    end
    %title(chName, 'FontSize',9, 'FontWeight','bold');
    pos_text = [xLims(1)+0.5*diff(xLims), yLims(1)+0.9*diff(yLims)];
    text(pos_text(1),pos_text(2), chName, ...
        'FontSize',9, 'FontWeight','bold', 'HorizontalAlignment','center');     % 
    
    % xlabel (bottom line)
    if ch > (nRows-1)*nCols
        if isfield(dataStruct.info, 'xlabel')
            xlabel(dataStruct.info.xlabel);
        end
    else
        set(gca, 'XTick',[]);
    end

    % ylabel (left column)
    if mod(ch,nCols) == 1
        if isfield(dataStruct.info, 'ylabel')
            ylabel(dataStruct.info.ylabel);
        end
    else
        set(gca, 'YTick',[]);
    end
    
    % vertical & horizontal bars (if specified)
    if isfield(dataStruct.info, 'verLines')
        for t = dataStruct.info.verLines
            plot([t, t], yLims, '--k');
        end
    end
    if isfield(dataStruct.info, 'horLines')
        for t = dataStruct.info.horLines
            plot([xLims(1), xLims(2)], [t, t], '--k');
        end
    end        
end

%% output directory
if isfield(dataStruct.info, 'outDir')
    outDir = dataStruct.info.outDir;
else
    outDir = params.storage.outputDir;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

%% figure name
if isfield(dataStruct.info, 'figName')
    figname = dataStruct.info.figName;
else
    figname = 'notNamed';
end

%% save
set(f, 'PaperPositionMode','auto');
saveas(f, [outDir filesep figname '.fig']);
print(f, '-dpng','-r0', [outDir filesep figname '.png']);
%print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname '.tif']);
close(f);    
display(['Figure: ' figname ' stored in: ' outDir]);

