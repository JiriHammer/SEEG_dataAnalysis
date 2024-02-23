function plotInfo = plotChannel_slicesActivity(params, D, H, selCh_H, plotInfo, groupInfo)
% plots into a single figure:
%   - brain slices (axial, sigittal,coronal) - also for bipolar montage
%   - freq. band activations
%   - spectra
% input: 
%   - D{freq} = trialsData saved from triggerAnalysis
%   - H = header
%   - selCh_H = selected channels (if bipolar = 2 x nCh)
%   - plotInfo = structure
% output:
%   - plotInfo: updated for fig, figName, outDir, was_significant (freq x 1)

% (c) Jiri, Jun20

%% time to plot
if ~isfield(groupInfo, 'time2plot')
    plotInfo.time2plot = [D{1}.xVals(1), D{1}.xVals(end)];  % all data
else
    plotInfo.time2plot = groupInfo.time2plot;       % specified data
end

%% variables
subjTag = plotInfo.subjTag;
freqBands = plotInfo.freqBands;
ch = plotInfo.ch;
gap = plotInfo.gap;
marg_h = plotInfo.marg_h;
marg_w = plotInfo.marg_h;
cacheFile = plotInfo.cacheFile;

%% figure: channel localization + freq. band activations
f = fig_make;
nSubs = 3*size(selCh_H,1) + size(D,1) + size(params.triggering.classes,1);
nRows = 3;
nCols = ceil(nSubs/nRows);
nSubs = nRows*nCols;
M_subs = reshape([1:nSubs], [nCols,nRows])';    % matrix of subplot numbers
v_subs = M_subs(:);                                  % vector of subplot numbers                        
plotInfo.fig = f;
plotInfo.nRows = nRows;
plotInfo.nCols = nCols;  
% plotInfo.marg_h = marg_h;
% plotInfo.marg_w = marg_h;
% plotInfo.gap = gap+0.03;
plotInfo.fontSize = 10;
plotInfo.labelSlices = false;

%% ---plot slices---
nPlot = 1;
for n = 1:size(selCh_H,1)                   % if n=2 -> bipolar
    thisCh = selCh_H(n,ch);                 % chnl index in H  
    plotInfo.thisCh = thisCh;
    plotInfo.channel = H.channels(thisCh);   

    % text on figure: chnl info
    plotInfo.txt_pos = [0.016, 0.98-0.03*(n-1)];
    textOnFig_chnlInfo(params, plotInfo);

    % subplot positions
    plotInfo.p_subPlot = v_subs([nPlot, nPlot+1, nPlot+2]);
    
    % slices (old)
%     plot_ieegChnlsInBrain(params, plotInfo);
    
    % slices (new)
    % set (or load) channel-specific brain volumes
    [brainVolumes, mni_vox] = brainVolumes_chSpec(params, params.plot_brainSlices, H.channels(thisCh), plotInfo.fixedVols, plotInfo.inds_vol, groupInfo);

    % >>> plot slices of selected channels (or MNI coors) <<<
    plotInfo.plot_brainSlices = params.plot_brainSlices;
    brainSlices_plot(brainVolumes, mni_vox, plotInfo);
    
    nPlot = nPlot+3;
end

%% ---plot activations---
was_significant = zeros(size(D,1),1);
for freq = 1:size(D,1)
    %subplot(nRows, nCols, nPlot);
    h_axes = subtightplot(nRows, nCols, v_subs(nPlot), 2*gap, 2*marg_h, 2*marg_w);
    hold on;
    plotInfo.chnlsAxes = h_axes;
    trialsData = D{freq};
    colors = trialsData.info.colors;
    trialsData.info.chNames{ch} = [];
    plotChannel2Axes(trialsData, ch, plotInfo);

    if v_subs(nPlot) > (nRows-1)*nCols
        xlabel('time [s]');
    end
    ylabel('rel.power');
    title(freqBands{freq,1});
    %-------------------------------------------------------------
    if isfield(trialsData, 'hVals')
%         if any(trialsData.hVals(:,ch) == 1)
%             was_significant(freq) = 1;
%         end
        was_significant(freq) = ch2roi_getSignificance(groupInfo, 'was_sgnf', trialsData.hVals(:,ch), trialsData.xVals);
    end
    nPlot = nPlot+1;
end
plotInfo.was_significant = was_significant;

%% ---spectra---
n = 1;
for clz = 1:size(params.triggering.classes,1)
    varName = ['spectralData_clz' num2str(clz)];
    % load spectra for given class
    if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
        disp(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
    else                
        clear(varName);
        load(cacheFile, varName);
        eval(['spectralData = ' varName ';']); 
        assert(size(selCh_H,2) == size(spectralData.cVals,3));

        % new: plot spectra
        plotInfo.cLims = [-0.5, 0.5];   % TO DO - automatic, old(now=May22) wo zscore: [-4, 4];
        plotInfo.nRows = nRows;
        plotInfo.nCols = nCols;
        plotInfo.marg_h = 2*marg_h;
        plotInfo.marg_w = 2*marg_w;
        plotInfo.gap = 2*gap;
        plotInfo.txt2axes = true;
        plotInfo.ylabel = 'freq (Hz)';
        spectralData.info.clr_clz = colors(clz,:); 
        plotInfo.thisSub = v_subs(nPlot);
        plotInfo.title_str = [params.triggering.classes{clz} '(' num2str(trialsData.info.nClz(n)) ')'];
        plotSpectrum2Axes(spectralData, ch, plotInfo);

        % colorbar
        if clz == size(params.triggering.classes,1)
            axPos = get(gca, 'Position');
            hCb = colorbar;
            ylabel(hCb, 'rel.power'); ... , 'FontSize',11, 'FontWeight','bold'
            set(hCb, 'position', [axPos(1)+1.1*axPos(3) axPos(2) 0.01 axPos(4)]);   
        end
        n = n+1;    % is this needed? n = clz ???
    end
    nPlot = nPlot+1;
end

%% save
figname=[subjTag '_' num2str(thisCh) '_' plotInfo.channel.name];
outDir =[params.storage.dir_results filesep params.storage.outName filesep 'slicesActivity_allChnls'];
fig_save(f, figname, outDir, 'format','png');      
plotInfo.figname = figname;
plotInfo.outDir = outDir;

