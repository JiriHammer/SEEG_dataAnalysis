function plotSpectrum2Axes(spectralData, ch, plotInfo)
% plots a spectrum of a channel 'ch' into axes

% (c) Jiri, Jun20

%% variables
cLims = plotInfo.cLims;
nRows = plotInfo.nRows;
nCols = plotInfo.nCols;
marg_h = plotInfo.marg_h;
marg_w = plotInfo.marg_w;
gap = plotInfo.gap;
thisSub =  plotInfo.thisSub;

if ~isfield(spectralData.info, 'clr_clz')
    clr = 'k';
else
    clr = spectralData.info.clr_clz;     % class-specific color
end

%% new subplot   
%subplot(nRows, nCols, ch);
subtightplot(nRows, nCols, thisSub, gap, marg_h, marg_w);
hold on;

%% >>> plot: imagesc <<<
h = imagesc(spectralData.xVals, spectralData.yVals, spectralData.cVals(:,:,ch), cLims);
% customColormap = brewermap(256,'*RdBu');
% colormap(customColormap);
colormap(gca,brewermap(256,'*RdBu'));   % set colormap for this axes (independent of figure colormap)

%% set y-limits
axis tight;
yLims = get(gca, 'ylim');
xLims = get(gca, 'xlim');

%% add also time-domain data
if isfield(spectralData, 'trial')
    set(gca, 'ylim', [-100, yLims(2)]);
    d = spectralData.trial.data(:,ch);
    d = linTransform(d, [min(d),max(d)], [-100, 0]);
    plot(spectralData.trial.time, d, 'k');
end

%% click on plot to see it bigger
set(gca, 'ButtonDownFcn','call_copy');
set(h, 'ButtonDownFcn','call_copy');

%% text: channel name
if isfield(plotInfo, 'txt2axes')
    if plotInfo.txt2axes
        if isfield(spectralData.info, 'chNames')
            chName = spectralData.info.chNames{ch};
        else
            chName = ['channel: ' num2str(ch)];
        end
        % break to 2 lines, if text is too long
        if size(chName,2) > 20
            i_semiColon = strfind(chName,':');
            if ~isempty(i_semiColon)
                tmp = chName;
                chName = cell(2,1);
                chName{1} = tmp(1:i_semiColon);
                chName{2} = tmp(i_semiColon+1:end);
            end
        end
        pos_text = [xLims(1)+0.5*diff(xLims), yLims(1)+0.9*diff(yLims)];
        text(pos_text(1),pos_text(2), chName, ...
            'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','center', 'Interpreter','none');     % 
    end
end

%% title
if isfield(plotInfo, 'title_str')
    title(plotInfo.title_str, 'Color',clr);   
end
%title([params.triggering.classes{clz} '(' num2str(trialsData.info.nClz(n)) ')'],'Color',colors(n,:));       % hot fix... (TO DO: trialsData.info.classes...) and remove 'n'    

%% xlabel (bottom line)
if thisSub > (nRows-1)*nCols
    if isfield(spectralData.info, 'xlabel')
        xlabel(spectralData.info.xlabel);
    end
else
    set(gca, 'XTick',[]);
end

%% ylabel (left column)
if isfield(plotInfo, 'ylabel')
    ylabel(spectralData.info.ylabel);
else
    if mod(ch,nCols) == 1
        if isfield(spectralData.info, 'ylabel')
            ylabel(spectralData.info.ylabel);
        end
    else
        set(gca, 'YTick',[]);
    end
end

%% vertical & horizontal bars (if specified)
plot([0, 0], yLims, ':k', 'LineWidth',2);
if isfield(spectralData.info, 'verLines')
    for t = spectralData.info.verLines
        if t > xLims(1) && t < xLims(2)
            plot([t, t], yLims, '--k');
        end
    end
end
if isfield(spectralData.info, 'horLines')
    for t = spectralData.info.horLines
        plot([xLims(1), xLims(2)], [t, t], '--k');
    end
end    

%% paradigm-specific times
if isfield(spectralData.info, 'paraTimes')
    if ~isempty(spectralData.info.paraTimes)
        for t = spectralData.info.paraTimes
            if ~isnan(t)
                if t > xLims(1) && t < xLims(2)
                    plot([t, t], yLims, 'Color',clr, 'LineWidth',2);
                end
            end
        end
    end
end

%% mark rejected indices
if isfield(spectralData.info, 'rejected')
    i_rej = find(spectralData.info.rejected(1,:,ch) == 1);
    for t = 1:length(i_rej)
        plot(spectralData.xVals(i_rej(t)), yLims(1)+0.7*diff(yLims), 'o', ...
            'MarkerSize',10, 'MarkerEdgeColor','k');    % , 'MarkerFaceColor','k'
    end        
end    
