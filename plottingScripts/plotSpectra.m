function plotSpectra(params, spectralData)
% plots spectra for channels
% spectralData = 
%      info: [1x1 struct]
%     xVals: [60x1 double]       ~ time
%     yVals: [129x1 double]      ~ freq
%     cVals: [129x60x125 double] ~ time x freq x chnls
    
% (c) Jiri, May16

%% basic check
nCh = size(spectralData.cVals,3);
assert(size(spectralData.cVals,1) == size(spectralData.yVals,1));
assert(size(spectralData.cVals,2) == size(spectralData.xVals,1));

%% color-limits
if isfield(spectralData.info, 'cLims')
    cLims = spectralData.info.cLims;
else
    cLims = getYLims(spectralData.cVals);
end

%% outliers
% if isfield(spectralData.info, 'outliers')
%     outliers = spectralData.info.outliers;
% else
%     outliers = getOutliers(spectralData.sVals);
% end
% assert(length(outliers) == nCh);

%% figure
f = fig_make;

%% figure layout
[nRows, nCols] = getSubplotLayout(nCh);
if ~isfield(spectralData.info, 'marg_h')
    marg_h = [0.05 0.05];
    marg_w = [0.04 0.04];
    gap = [0.005, 0.005];
else
    marg_h = spectralData.info.marg_h;
    marg_w = spectralData.info.marg_w;
    gap = spectralData.info.gap;
end

%% text on upper part of the figure
if isfield(spectralData.info, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = spectralData.info.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold');
end

%% settings
plotInfo.nCh = nCh;
plotInfo.cLims = cLims;
plotInfo.nRows = nRows;
plotInfo.nCols = nCols;
plotInfo.marg_h = marg_h;
plotInfo.marg_w = marg_w;
plotInfo.gap = gap;
plotInfo.txt2axes = true;

%% >>> plot channel-wise <<<
for ch = 1:nCh
    plotInfo.thisSub = ch;
    plotSpectrum2Axes(spectralData, ch, plotInfo);
end

%% colorbar
if isfield(spectralData.info, 'clrbar_title')
    clrbar_title = spectralData.info.clrbar_title;
else
    clrbar_title = 'rel.power (dB)';
end
hCb = colorbar;
ylabel(hCb, clrbar_title, 'FontSize',11, 'FontWeight','bold');
set(hCb, 'position', [0.97 0.2 0.01 0.6]);

%% figure properties
set(gcf,'PaperPositionMode','auto');
set(gcf, 'color', [1 1 1]);   %white background

%% output directory
if isfield(spectralData.info, 'outDir')
    outDir = spectralData.info.outDir;
else
    outDir = params.storage.outputDir;
end

%% figure name
if isfield(spectralData.info, 'figName')
    figname = spectralData.info.figName;
else
    figname = 'notNamed';
end

%% save
fig_save(f, figname, outDir);
close(f);    

