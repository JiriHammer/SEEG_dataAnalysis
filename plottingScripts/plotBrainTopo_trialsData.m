function plotBrainTopo_trialsData(params, trialsData)
% plots trialsData onto brain topology (slices or 3D model)
% trialsData = 
%      info: [1x1 struct]
%     xVals: [1025x1 double]        ~ 1D: time, in [s]
%     yVals: [1025x110x2 double]    ~ 3D: time x channels x classes
%     yErrs: [1025x110x2 double]    ~ 3D: time x channels x classes

% (c) Jiri, May17

%% default settings
if ~isfield(trialsData.info, 'name_freqBand')
    name_freqBand = 'noFreqBandName';
else
    name_freqBand = trialsData.info.name_freqBand;
end

if ~isfield(trialsData.info, 'mriFile_slices')
    mriFile_slices = 'wT1';
else
    mriFile_slices = trialsData.info.mriFile_slices;
end

if ~isfield(trialsData.info, 'mriFile_3Dmodel')
    mriFile_3Dmodel = 'wc1T1';
else
    mriFile_3Dmodel = trialsData.info.mriFile_3Dmodel;
end

if ~isfield(trialsData.info, 'chnlsMNI')
    chnlsMNI = [];
else
    chnlsMNI = trialsData.info.chnlsMNI;
end

if ~isfield(trialsData.info, 'outDir_beg')
    outDir_beg = params.storage.outputDir;
else
    outDir_beg = trialsData.info.outDir_beg;
end

%% plot brain topo in slices
if params.plot_triggering.brainTopo && params.plot_brainTopo.plot_slices
    % brain info
    plotInfo.plottingStyle = 'slices';
    plotInfo.file2load = mriFile_slices;      % choices: wrCT, wT1, wT2, wc1T1, wc2T2, ...
    plotInfo.brain = getBrainData(params, plotInfo);
    
    % channels MNI coors
    if isempty(chnlsMNI)
        plotInfo.chnlsMNI = getSelChnlsMNI(params);
    else
        plotInfo.chnlsMNI = chnlsMNI;
    end
    
    % plot
    time2plot = params.plot_brainTopo.time2plot;        % in [s], w.r.t. cutting point 
    clims = [prctile(trialsData.yVals(:),1), prctile(trialsData.yVals(:),99)]; clims = max(abs(clims)); 
    plotInfo.chnl_clims = [-clims, clims];
    for clz = 1:size(trialsData.yVals,3)
        for t = 1:size(time2plot,2)
            i_t = closestval(trialsData.xVals, time2plot(t));
            vals = trialsData.yVals(i_t,:,clz);
            plotInfo.text = ['subject: ' params.storage.subjTag ...
                             ', triggered: ' params.triggering.cutPoint ...
                             ', class: ' trialsData.info.clzNames{clz,1} '(' num2str(trialsData.info.nClz(clz)) ')' ...
                             ', freq. band: ' name_freqBand ...
                             ', time = ' num2str(time2plot(t)) ' s'];   
            plotInfo.outDir =[outDir_beg filesep 'topoBrain_' name_freqBand filesep 'class_' num2str(clz) '_' trialsData.info.clzNames{clz,1}];
            plotInfo.figName=['fig_' num2str(i_t)];
            plot_brainSlices(params, vals(:), plotInfo);
        end
    end
end

%% plot brain topo in 3D model
if params.plot_triggering.brainTopo && (params.plot_brainTopo.plot_projections || params.plot_brainTopo.plot_animation)
    % brain info
    plotInfo.plottingStyle = '3D_model';
    plotInfo.file2load = mriFile_3Dmodel;      % choices: wrCT, wT1, wT2, wc1T1, wc2T2, ...
    plotInfo.brain = getBrainData(params, plotInfo);
    
    % anatomy areas
    if isfield(trialsData.info, 'mriFile_AnatAreas')
        plotInfo.aarea = cell(size(trialsData.info.mriFile_AnatAreas,1),1);
        for aa = 1:size(trialsData.info.mriFile_AnatAreas,1)
            plotInfo.file2load = trialsData.info.mriFile_AnatAreas{aa,1};      % choices: see: TO DO
            plotInfo.aarea{aa} = getBrainData(params, plotInfo);
            plotInfo.aarea{aa}.faceClr = trialsData.info.mriFile_AnatAreas{aa,2};
        end
    end
    
    % channels MNI coors
    if isempty(chnlsMNI)
        plotInfo.chnlsMNI = getSelChnlsMNI(params);
    else
        plotInfo.chnlsMNI = chnlsMNI;
    end
    
    % plot
    time2plot = params.plot_brainTopo.time2plot;        % in [s], w.r.t. cutting point 
    plotInfo.clims = [prctile(trialsData.yVals(:),1), prctile(trialsData.yVals(:),99)];
    plotInfo.fig_bkgClr = 'k';                          % choices: black ('k') or white ('w')
    for clz = 1:size(trialsData.yVals,3)
        for t = 1:size(time2plot,2)
            i_t = closestval(trialsData.xVals, time2plot(t));
            vals = trialsData.yVals(i_t,:,clz);
            plotInfo.text = ['subject: ' params.storage.subjTag ...
                             ', triggered: ' params.triggering.cutPoint ...
                             ', class: ' trialsData.info.clzNames{clz,1} '(' num2str(trialsData.info.nClz(clz)) ')' ...
                             ', time = ' num2str(time2plot(t)) ' s'];   
            plotInfo.outDir =[outDir_beg filesep 'topoBrain_' name_freqBand filesep 'class_' num2str(clz) '_' trialsData.info.clzNames{clz,1}];
            plotInfo.figName=['fig_' num2str(i_t)];
            if params.plot_brainTopo.plot_animation
                plot_brain3D_wrapper(params, vals(:), plotInfo);
            end
            if params.plot_brainTopo.plot_projections
                brain3D_projections(params, vals(:), plotInfo);
            end
        end
    end
end
