function plotHist(params, histData)
% plots histogram distribution of 'histData'
% histData = 2D: samples x chnls (may inlcude NaNs)

% (c) Jiri, Aug16

%% load 'H'
load(params.storage.cacheFile, 'H');
nCh = size(histData,2);
nBins = round(size(histData,1)/100);

% channel names
params.nCh = size(histData,2);
chNames = getChannelNames(params, params.response.signalType);
    
%% figure
f = figure('units','normalized','outerposition',[0 0 1 1]);

[nRows, nCols] = getSubplotLayout(nCh);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.020, 0.010];

%% text on upper part of the figure
if isfield(params, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = params.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fonts', 16, 'fontw', 'bold');
end

%% plot: channels distrubution
for ch = 1:size(histData,2)  
    subtightplot(nRows, nCols, ch, gap, marg_h, marg_w);
    hold on;
    hist(histData(:,ch,:), nBins);
%     h = findobj(gca,'Type','patch');
%     set(h,'EdgeColor','none');

    % set wider x-limits
    axis tight;
    xLims = get(gca, 'xlim');
    set(gca, 'xlim',[xLims(1)-0.05*diff(xLims), xLims(2)+0.05*diff(xLims)]);
    yLims = get(gca, 'ylim');
    
    % title = channel name
    %title(chNames{ch}, 'Interpreter','none');
    pos_text = [xLims(1)+0.5*diff(xLims), yLims(1)+0.9*diff(yLims)];
    chName = chNames{ch};
    if size(chName,2) > 20
        i_semiColon = strfind(chName,':');
        if ~isempty(i_semiColon)
            tmp = chName;
            chName = cell(2,1);
            chName{1} = tmp(1:i_semiColon);
            chName{2} = tmp(i_semiColon+1:end);
        end
    end    
    text(pos_text(1),pos_text(2), chName, ...
        'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','center', 'Color','k', 'Interpreter','none');     % 
        
    % plot mean & median
    m = nanmean(histData(:,ch),1);
    plot([m, m], yLims, '--r', 'LineWidth',2);
    m = nanmedian(histData(:,ch),1);
    plot([m, m], yLims, '--m');
    
    % vertical & horizontal bars (if specified)
    if isfield(params, 'verLines')
        for t = params.verLines
            plot([t, t], yLims, '--k');
        end
    end
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn','call_copy');    
end

%% save
if isfield(params, 'outDir')
    outDir = params.outDir;
else
    outDir = [params.storage.outputDir filesep 'tuning_respDistr'];
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

if isfield(params, 'figName')
    figname = params.figName;
else
    figname = 'resp_sessDistr';
end

set(f, 'PaperPositionMode','auto');
saveas(f, [outDir filesep figname '.fig']);
print('-dpng','-r0', [outDir filesep figname '.png']);
disp(['Figure: ' figname ' stored in: ' outDir]);
close(f);   

