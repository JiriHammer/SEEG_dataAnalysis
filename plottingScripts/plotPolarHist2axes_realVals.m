function plotPolarHist2axes_realVals(ampPhiBins, plotInfo)
% plots the tuning bins in: ampPhiBins [ampBin x phiBin] of real values
% into polar histogram, axes position spec. in 'plotInfo'

% (c) Jiri, Mar12

%% plot colorbar
if plotInfo.colorbar
    % vals
    nCols = 2;
    nRows = 359;
    ...aRes = 1/nCols;
    aVector = 1:nCols;
    cRes = sum(abs(plotInfo.clims))/nRows;
    cVector = [plotInfo.clims(1):cRes:plotInfo.clims(2)]';
    ...alphaVals = repmat(aVector, size(cVector,1), 1);
    colorVals = repmat(cVector,1,size(aVector,2));
    X = aVector;
    Y = cVector;
    Z = ones(size(colorVals,1),size(colorVals,2));

    % axes
    ax = axes('Parent',plotInfo.figHandle, 'Position',plotInfo.axPos, 'NextPlot','replacechildren');
    h = pcolor(X,Y,Z);
    axis tight 
    set(h, ...
...        'AlphaDataMapping','none', ...
...        'FaceAlpha','interp', ...
...        'EdgeAlpha','interp', ...    
...        'AlphaData',alphaVals,... 
        'EdgeColor','none', ...        %'interp', ...
        'CData',colorVals ...
        );
    if isfield(plotInfo, 'cVals2cInds')         % for diff. colormaps in 1 figure
        cInds = cVals2cInds(colorVals, [-pi, pi], plotInfo.mapInds);
        set(h, 'CData',round(cInds));
        caxis([1 size(get(gcf,'ColorMap'),1)]);
    else
        colormap(ax, plotInfo.clrMap);
        set(gca, 'CLim',plotInfo.clims);
    end
    set(gca,'ButtonDownFcn','call_copy'); %click on plot to see it bigger
    
    % annotation
    if isfield(plotInfo, 'ytickVals'), set(gca, 'YTick',plotInfo.ytickVals); end
    if isfield(plotInfo, 'ytickLabels'), set(gca,'YTickLabel',plotInfo.ytickLabels); end
    if isfield(plotInfo, 'yLabel'), set(gca,'YLabel',plotInfo.yLabel); end
    set(ax, 'XTick',[]);
%     xTickVals = [aVector(1), aVector(round(length(aVector)/2)), aVector(end)];
%     set(gca,'XTick',xTickVals);
%     xTickLabels = cell(1,length(xTickVals));
%     for label = 1:length(xTickVals)
%         xTickLabels{label} = num2str(xTickVals(label), '%01.2f');
%     end
%     set(gca,'XTickLabel',xTickLabels);
%    xlabel('amplitude');
%    set(gca, 'Box','on', 'LineWidth',1.5);
    return;
end

%% plot nEvents (but also auxChannels, etc.)
if plotInfo.nEvents
    % axes
    ax = axes('Parent',plotInfo.figHandle, 'Position',plotInfo.axPos, 'NextPlot','replacechildren');
    axis square;
    set(ax, 'XTick',[], 'YTick',[]);
    set(ax, 'XColor','k', 'YColor','k', 'Box','on', 'LineWidth',1.5);
    hold on;

    % increase size of ampPhiBins (pcolor ignores last row and last col)
    ampPhiBins = cat(2, ampPhiBins, ampPhiBins(:,end));     % duplicate last col 
    ampPhiBins = cat(1, ampPhiBins, ampPhiBins(end,:));     % duplicate last row 
    
    % bin vals
    %ampR = plotInfo.ampBinVals';    % bin values
    ampR = [1:length(plotInfo.ampBinVals)]';      % bin index
    phiA = plotInfo.phiBinVals;
    assert(size(ampR,1) == size(ampPhiBins,1));
    assert(size(phiA,2) == size(ampPhiBins,2));

    X = ampR*cos(phiA);
    Y = ampR*sin(phiA);

    % transparency & color data
    alphaData = ones(size(X));
    colorData = ampPhiBins;
    Z = ones(size(X));

    % pcolor - plot
    h = pcolor(X,Y,Z);
    axis equal tight 
    %colormap(ax, 'hsv');
    set(h, ...
...        'AlphaDataMapping','none', ....
...        'FaceAlpha','interp', ...
...        'EdgeAlpha','interp', ...    
...        'AlphaData',alphaData,... 
        'EdgeColor','none', ...        %'interp', ...
        'CData',colorData ...
        );
    set(gca,'xtick',[],'ytick',[],'CLim',plotInfo.clim);
    set(h,'ButtonDownFcn','call_copy');                     % click on plot to see it bigger
    set(gca,'ButtonDownFcn','call_copy');                   % click on plot to see it bigger
    title(plotInfo.title, 'FontSize',12, 'fontw','bold');
    colorbar;
    return;
end

%% channels: axes position
ax = axes('Parent',plotInfo.figHandle, 'Position',plotInfo.axPos, 'NextPlot','replacechildren');
axis square;
set(ax, 'XTick',[], 'YTick',[]);
set(ax, 'XColor',plotInfo.aaClr, 'YColor',plotInfo.aaClr, 'Box','on', 'LineWidth',1.5);
hold on;

%% increase size of ampPhiBins (pcolor ignores last row and last col)
ampPhiBins = cat(2, ampPhiBins, ampPhiBins(:,end));     % duplicate last col 
ampPhiBins = cat(1, ampPhiBins, ampPhiBins(end,:));     % duplicate last row 

%% bin vals
%ampR = plotInfo.ampBinVals';    % bin values
ampR = [1:length(plotInfo.ampBinVals)]';      % bin index
phiA = plotInfo.phiBinVals;
assert(size(ampR,1) == size(ampPhiBins,1));
assert(size(phiA,2) == size(ampPhiBins,2));

%% channels: vals & transparency & color data
X = ampR*cos(phiA);
Y = ampR*sin(phiA);
Z = ones(size(X));

%alphaData = transparencyCoding(ampCh, 'code', plotInfo); %ampCh./plotInfo.prctiles;                           % linear mapping, transparency values: [0, 1]
colorData = ampPhiBins;

%% channels: pcolor - plot
h = pcolor(X,Y,Z);
axis equal tight 
set(h, ...
...    'AlphaDataMapping','none', ...
...    'FaceAlpha','interp', ...       % choice: 'flat', 'interp'
...    'EdgeAlpha','interp', ...    
...    'AlphaData',alphaData,... 
    'EdgeColor','none', ...        %'interp', ...
    'CData',colorData ...
    );
%shading(gca,'interp');
set(gca,'xtick',[],'ytick',[]);
set(h,'ButtonDownFcn','call_copy');                     % click on plot to see it bigger
set(gca,'ButtonDownFcn','call_copy');                   % click on plot to see it bigger
    
if isfield(plotInfo, 'cVals2cInds')         % for diff. colormaps in 1 figure
    cInds = cVals2cInds(colorData, plotInfo.clims, plotInfo.mapInds);
    set(h, 'CData',cInds);
    caxis([1 size(get(gcf,'ColorMap'),1)]);
else
    colormap(ax, plotInfo.clrMap);
    set(gca, 'CLim',plotInfo.clims);
end

% title
if isfield(plotInfo, 'title')
    title(plotInfo.title, 'FontSize',12, 'fontw','bold');
end
