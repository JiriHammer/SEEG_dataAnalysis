function topoplotMeanWithError(params, statData)
% plots as a grid values in 'statData' cell struct with fields:
%    yVals: [sample x ch]
%    yErrs: [sample x ch]
%    xVals: x-Axis values (vector)
%    info: struct with supplementary informatio

% (c) Jiri, Jan11

for sig = 1:size(statData,1)
    assert(size(statData{sig,1}.yVals,1) == length(statData{1,1}.xVals));
    
    % channel names and pos
    signalType = getStructValForFieldname(statData{sig,1}.info, 'signalType');
    chNamesAndPos = getChNamesAndPos(params.connectionTable, signalType);
    chNames = getChannelNames(params.connectionTable, signalType);
    selCh = getStructValForFieldname(statData{1,1}.info, 'selCh');
    chNames = chNames(selCh);
    
    % auxiliary information of plot
    titleSuffix = getStructValForFieldname(statData{1,1}.info, 'title');
    
    % y-lim
    if isfield(statData{1,1}.info, 'yLims')
        ylimits = statData{1,1}.info.yLims;
        outliers = [];
    else
        if isfield(statData{1,1}, 'yErrs')
            inclErrors = true;
        else
            inclErrors = false;
        end
        [ylimits, outliers] = getYLims(statData(sig,:), inclErrors);
    end
    
    % figure
    f = figure('visible', 'on', 'Position', [1, 1, 1920, 1200]);
    set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);

    % text
    tx =axes('visible','off', 'position',[0 0 1 1]);
    if isfield(statData{sig,1}.info, 'text')
        mytitle = statData{sig,1}.info.text;
    else
        mytitle = ['mean +/- std over sessions. signalType: ' signalType];
    end
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fonts', 14, 'fontw', 'bold'); 
    
    % plot channel-wise
    assert(size(statData{sig,1}.yVals,2) == length(chNames));
    for ch = 1:size(statData{sig,1}.yVals,2)
        % channel name and position
        chName = chNames{ch};
        if strcmp(signalType, 'EcogGrid') && ~isempty(findstr(params.job.experimentId, 'WN'))
            aaClr = getChClrOfAa(params.connectionTable, signalType, selCh(ch));
        elseif strcmp(signalType, 'EcogStrip') && ~isempty(findstr(params.job.experimentId, 'CS'))
            aaClr = getClrForAa(chName);
        else
            aaClr = [0 0 0];  % ~ black
        end
        [left bottom width height] = getBoxPosForName(chNamesAndPos, chName);
        
        % detect outlier and set y-limits
        if ismember(ch, outliers)
            bkClr = 'r';
        else
            bkClr = 'w';
        end
        
        % new subplot
        axes('OuterPosition',[left bottom width height], 'visible','on', 'Color',bkClr);
        set(gca, 'XColor',aaClr, 'YColor',aaClr, 'Box','on', 'LineWidth',1.5);        
        hold on;
        
        % mean - error band plot
        for clz = 1:size(statData,2)
            if isfield(statData{1,1}, 'yErrs')
                h(:,clz) = plotband(statData{sig,clz}.xVals, statData{sig,clz}.yVals(:,ch), statData{sig,clz}.yErrs(:,ch), statData{sig,clz}.info.clr);
            else
                h(:,clz) = plot(statData{sig,clz}.xVals, statData{sig,clz}.yVals(:,ch), 'Color',statData{sig,clz}.info.clr, 'LineWidth',1);
            end
        end
        axis tight;
        
        % set y-limits
        if ismember(ch, outliers)
            chYLimits = get(gca, 'yLim');
        else
            chYLimits = ylimits;
        end  
        if ~ischar(ylimits)     % hack: use character to suppress ylims computation...
            ylim(chYLimits);
        else
            chYLimits = get(gca, 'yLim');
            set(gca, 'ylim', [chYLimits(1)-0.2*diff(chYLimits), chYLimits(2)+0.2*diff(chYLimits)]);
            chYLimits = get(gca, 'yLim');
        end
        title([chName '(' num2str(selCh(ch)) ')' titleSuffix]);
        if ch == 1
            xlabel(statData{1,1}.info.xlabel);
        end
        
        % vertical & horizontal bars (if specified)
        if isfield(statData{1,1}.info, 'verLines')
            for t = statData{1,1}.info.verLines
                plot([t, t], chYLimits, '--k');
            end
        end
        if isfield(statData{1,1}.info, 'horLines')
            xLims = get(gca, 'xlim');
            for t = statData{1,1}.info.horLines
                plot([xLims(1), xLims(2)], [t, t], '--k');
            end
        end        
    end

    % output directory
    if isfield(statData{sig,1}.info, 'outDir')
        outDir = statData{sig,1}.info.outDir;
    else
        outDir = params.storage.outputDir;
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end  

    % save
    figname = getStructValForFieldname(statData{sig,1}.info, 'figname');
    saveas(f, [outDir filesep figname '.fig']);
    print(f, '-dpng','-r0', [outDir filesep figname '.png']);
    %print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname]);
    %close(f);    
    display(['Figure: ' figname ' stored in: ' outDir]);
end