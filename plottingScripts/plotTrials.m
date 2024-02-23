function plotTrials(params, trialsData)
% plots values in 'trialsData' struct with fields:
%    yVals: [sample x ch x clz]
%    yErrs: [sample x ch x clz] -> optional (plotband)
%    xVals: x-Axis values (vector), typically time in [s]
%    info: struct with supplementary information

% (c) Jiri, Jan11 (in whitenoise), Mar16

%% basic check
nCh = size(trialsData.yVals,2);
assert(size(trialsData.yVals,1) == size(trialsData.xVals,1));

%% settings
if isfield(trialsData.info, 'fontsize')
    fontsize = trialsData.info.fontsize;
else
    fontsize = 12;
end

%% load 'H'
%load(params.storage.cacheFile, 'H');

%% y-lim
if isfield(trialsData.info, 'yLims')
    yLims = trialsData.info.yLims;      % if empty, no common y-limits!
else
    if isfield(trialsData, 'yErrs')
        yLims_lo = getYLims(trialsData.yVals - trialsData.yErrs);
        yLims_hi = getYLims(trialsData.yVals + trialsData.yErrs);
        yLims = [yLims_lo(1), yLims_hi(2)];
    else
        yLims = getYLims(trialsData.yVals);
    end
end

%% outliers
if isfield(trialsData.info, 'outliers')
    outliers = trialsData.info.outliers;
else
    outliers = getOutliers(trialsData.yVals);
end
assert(length(outliers) == nCh);
%outliers = zeros(nCh,1);    % hack!

%% paradigm times: vertical lines (reaction times)
% if isfield(params.plot_triggering, 'paraTimes')
%     if ismember('RT', params.plot_triggering.paraTimes)
%         assert(ismember('RT', params.triggering.behaviour));    % 'RT' must be precomputed & saved
%         
%         % load reaction times 'RT' from cacheFile
%         clear RT
%         load(params.storage.cacheFile, 'RT');
%         assert(size(RT,1) == size(trialsData.yVals,3));
%         trialsData.info.RT = RT;
%     end
% end

%% figure
f = fig_make;

[nRows, nCols] = getSubplotLayout(nCh);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.006, 0.006];

%% text on upper part of the figure
if isfield(trialsData.info, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = trialsData.info.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 16, 'fontw', 'bold');
end

%% plot channel-wise
for ch = 1:nCh

    % if outlier, set background color to pink
    if any(outliers(ch,:))
        bkgClr = [1 0.8 0.8];
        txtClr = 'r';
    else
        bkgClr = [1 1 1];
        txtClr = 'k';
    end

    % new subplot   
    %h = subplot(nRows, nCols, ch);
    h = subtightplot(nRows, nCols, ch, gap, marg_h, marg_w);
    set(gca, 'Color',bkgClr);
    hold on;

    % >>>>>>>>>>>>>>> plot <<<<<<<<<<<<<<<<<<<<
    plotInfo.chnlsAxes = h;
    plotInfo.yLims = yLims;
    plotInfo.txtClr = txtClr;
    wasOutOfScale = plotChannel2Axes(trialsData, ch, plotInfo);
%     drawnow;
    % <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>
    
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
            ylabel(trialsData.info.ylabel, 'FontSize',fontsize);
        end
    else
        if wasOutOfScale  % set y-scale
            1+1;
        elseif ~isempty(yLims)  % || ~outliers(ch,:)
            set(gca, 'YTick',[]);   % do not set y-scale
        end
    end
    set(gca, 'FontSize',fontsize);
    
    % uncomment for debugging
%     disp(['done plotting ch = ' num2str(ch)]);
%     if ch == 186 || ch == 1
%         why;
%     end
end

%% output directory
if isfield(trialsData.info, 'outDir')
    outDir = trialsData.info.outDir;
else
    outDir = params.storage.outputDir;
end

%% figure name
if isfield(trialsData.info, 'figName')
    figname = trialsData.info.figName;
else
    figname = 'notNamed';
end

%% save
fig_save(f, figname, outDir, 'format','png');
close(f);    

