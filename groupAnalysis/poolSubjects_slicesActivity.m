function poolSubjects_slicesActivity(params, groupInfo)
% plots brain slices (ax,cor,sag) & activations (ERP + freq.band + spectra)
% pools together subject processed by a given job
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Apr18

if params.plot_brainSlices.plot_slices
    disp('Plotting channels: 3D slices & activations ...');
else
    disp('Plotting channels: 3D slices & activations: NOT selected.');
    return;
end

%% settings: default
if nargin < 2
    groupInfo = struct;
end

%% settings: activations
freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands

%% load & plot frequency bands + ERP (over subjects and anatomy areas)
plotInfo.freqBands = freqBands;
plotInfo.subtractAvgResponse = groupInfo.subtractAvgResponse;
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    params.storage.subjTag = subjTag;
    plotInfo.subjTag = subjTag;
    
    %% cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    plotInfo.cacheFile = cacheFile;

    %% load freq. band activations
    D = cell(size(freqBands,1),1);
    for freq = 1:size(freqBands,1)
        varName = ['trialsData_' freqBands{freq,1}];  
        % try to load freq. band activation from cache file
        if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
            error(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
        else
            load(cacheFile, varName);
            assert(exist(varName,'var') == 1);
            eval(['trialsData = ' varName ';']); 
            clear(varName); 
            D{freq} = trialsData;
            clear trialsData;
        end
    end
  
    %% load structure H & selected channels
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');    
    selCh_H = selCh_H_resp;

    %% anatomical assignements from isarg_atlas
    [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, cacheFile);

    %% load brain data
    [fixedVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brainSlices);
    plotInfo.clrmap = clrmap;   
    plotInfo.fixedVols = fixedVols;
    plotInfo.inds_vol = inds_vol;

    %% plot channel activations & slices
    for ch = 1:size(selCh_H,2)
        
        %% >>> plot channel: slices & activations (saves fig to: slicesActivity) <<<
        plotInfo.ch = ch;
        plotInfo.marg_h = [0.05 0.05];
        plotInfo.marg_w = [0.04 0.04];
        plotInfo.gap = [0.03, 0.03];
        plotInfo = plotChannel_slicesActivity(params, D, H, selCh_H, plotInfo, groupInfo);
        f = plotInfo.fig;

        %% save figure to different directories (as significant, ...)
        groupInfo.path_isarg_atlas = params.path2others.isarg_atlas;
        anatomy_signficance_saveFig(params, plotInfo, groupInfo, ass_isargAtlas(ch));

        %% close fig
        close(f);    
    end
end

