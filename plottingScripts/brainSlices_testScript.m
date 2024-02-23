%% plot slices: test script 
%   - load header H
%   - select iEEG channels
%   - load brain volumes that are "fixed" (same for all channels of one subject, e.g. wT1.nii)
%       -> brainVolumes_fixed.m
%   - load brain volumes that are "channel-specific" (e.g., wDefault.nii)
%       -> brainVolumes_chSpec.m
%   - plot slices
%       -> brainSlices_plot.m

% (c) Jiri, Nov21

%% settings
params_default;
params.storage.subjTag = '20_PR4';
params.storage.pathBeg = 'G:\dox\ms_switch_EI\data\v2_stft_sessions\switchin_EI_IE_car';

params.plot_brainSlices = struct(...           % user interface structure: holds most (but not all!) of the user settings
    'plot_slices', true, ...                            % 2D brain slices
    'volumes2plot', {{...       % volumes to plot (coregistered .nii files)
        'wT1', gray(128), [1, 1];  ...     % name & colormap & transparency limits
        'wrCT', copper(128), [0 0.8]; ...
        'ch-spec_Yeo7', autumn(128), [0, 0.4]; ...
        }}, ...             
    'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'size_interpolate', 1.0, ...                        % in [mm], voxel size to which the brain is interpolated
    'printFormats', {{'fig', 'png'}}, ...
    'printResolution', 0 ...                           % choices: 0 (= screen resolution) or 600 (= dpi). Resolution of the figures.                          
); 

%% ============================ IMPLEMENTATION ============================

%% load header H
fileName = [params.storage.pathBeg filesep params.storage.subjTag filesep 'cacheFile.mat'];
load(fileName, 'H');

%% selected iEEG channels
i_sel = [];
for ch = 4:size(H.channels,2)
    if ismember(H.channels(ch).signalType, {'SEEG', 'ECoG-Grid', 'ECoG-Strip'})
        i_sel = cat(2, i_sel, ch);
    end
end
channels = H.channels(i_sel);

%% load brain data -> fixedVols & clrmap
[fixedVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brainSlices);
plotInfo.clrmap = clrmap;

%% update ch-specific brain volumes & plot slices at channel's MNI coors 
plotInfo.outDir = ['G:\dox\ms_switch_EI\data\v2_stft_sessions\switchin_EI_IE_car\test_slices' filesep 'plot_channelLocalizations'];
plotInfo.subjTag = params.storage.subjTag;
for ch = 1:size(channels,2)
    
    % set (or update) channel-specific brain volumes
    [brainVolumes, mni_vox] = brainVolumes_chSpec(params, params.plot_brainSlices, channels(ch), fixedVols, inds_vol);

    % >>> plot slices of selected channels (or MNI coors) <<<
    plotInfo.channel = channels(ch);
    plotInfo.thisCh = i_sel(ch);        % channel index in H.channels
    plotInfo.plot_brainSlices = params.plot_brainSlices;
    brainSlices_plot(brainVolumes, mni_vox, plotInfo);
end

