function plotChannel2Axes_complexBins(ampPhiBins, plotInfo)
% plots the tuning bins in: ampPhiBins [ampBin x phiBin]
% uses polar coordinates (amp x phi bins) for each channel
% complex predictor (e.g. 2D velocity), real response (e.g. high-gamma iEEG channel)
% axes position defined in 'plotInfo'

% ampPhiBins = 2D: [amp x phi] bins
% plotInfo = 
%            ampBinVals: [0 5 10 20 30 60 150]
%         ampBinCenters: [6×1 double]
%            phiBinVals: [-3.1416 -2.3562 -1.5708 -0.7854 0 0.7854 1.5708 2.3562 3.1416]
%         phiBinCenters: [8×1 double]
%               thisLag: -2.4000
%               figname: 'lag_3'
%     plotSmallColorbar: 0
%             figHandle: [1×1 Figure]
%                 axPos: [0.0400 0.8257 0.0978 0.1243]
%                 aaClr: 'k'
%                 cLims: [-0.7274 0.6949]
%              colorbar: 0
             
% (c) Jiri, Sep18
% based on: whitenoise/plotting/plotChBins2circularAxes.m (Feb12)

%% default settings
if ~isfield(plotInfo, 'colorbar'), plotInfo.colorbar = false; end
if ~isfield(plotInfo, 'nEvents'), plotInfo.nEvents = false; end
if ~isfield(plotInfo, 'plotSmallColorbar'), plotInfo.plotSmallColorbar = false; end

%% plot big colorbar
if plotInfo.colorbar
    axClrbar = colorbar('location','EastOutside', 'position', [0.97 0.2 0.01 0.6], 'FontSize',10);
    set(get(axClrbar,'YLabel'),'String','resp [std]', 'FontSize',10);
    return;
end

%% axes position
ax = axes('Parent',plotInfo.figHandle, 'Position',plotInfo.axPos, 'NextPlot','replacechildren');
axis square;
set(ax, 'XTick',[], 'YTick',[]);
set(ax, 'XColor',plotInfo.aaClr, 'YColor',plotInfo.aaClr, 'Box','on', 'LineWidth',0.5);
% if isfield(plotInfo, 'axBkgClr')
%     set(ax, 'Color', plotInfo.axBkgClr);
% end
hold on;

%% axes vals
phiA = plotInfo.phiBinVals;
ampR = [0:size(plotInfo.ampBinVals,2)-1]';  ... plotInfo.ampBinVals'
X = ampR*cos(phiA);
Y = ampR*sin(phiA);
C = ampPhiBins;
C = cat(2, C, C(:,end));    % add last col (dummy variable)
C = cat(1, C, C(end,:));    % add last row (dummy variable)
Z = ones(size(X));

%% polar plot: pcolor
h = pcolor(X,Y,Z);
axis equal tight 
set(h, ...
    'EdgeColor','none', ...        % opt: 'interp', ...
    'CData',C ...
    );
shading flat;
set(gca,'xtick',[],'ytick',[]);

%% colormaps
if isfield(plotInfo, 'cVals2cInds')         % for diff. colormaps in 1 figure
    cInds = cVals2cInds(C, [-pi, pi], plotInfo.mapInds);
    set(h, 'CData',cInds);
    caxis([1 size(get(gcf,'ColorMap'),1)]);
else
    colormap(ax, 'jet');
    %colormap(ax, getColorMap('bwwwr'));
    set(gca, 'CLim',plotInfo.cLims);
end

%% axes has its own little colorbar?
if plotInfo.plotSmallColorbar
    colorbar;
end

