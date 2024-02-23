function brainSlices_prep2plot(params, channels, plotInfo)
% brain slices plotting: prepares necessary data
% channels = H.channels(i_sel)

% (c) Jiri, Nov21
% based on plot_ieegChnlsInBrain.m

%% plot brain slices?
if ~isfield(params, 'plot_brainSlices')
    params_default;
end
if ~params.plot_brainSlices.plot_slices
    disp('plotting of brain sluices not selected.');
    return;
end

%% load brain data -> fixedVols & clrmap
[fixedVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brainSlices);
plotInfo.clrmap = clrmap;

%% plot slices at channel's MNI coors 
for ch = 1:size(channels,2)
    
    % set (or load) channel-specific brain volumes
    [brainVolumes, mni_vox] = brainVolumes_chSpec(params, params.plot_brainSlices, channels(ch), fixedVols, inds_vol);

    % >>> plot slices of selected channels (or MNI coors) <<<
    plotInfo.channel = channels(ch);
    plotInfo.thisCh = plotInfo.selChannels(ch);
    plotInfo.plot_brainSlices = params.plot_brainSlices;
    plotInfo.figName=[num2str(ch) '_' plotInfo.channel.name];
    brainSlices_plot(brainVolumes, mni_vox, plotInfo);
end

