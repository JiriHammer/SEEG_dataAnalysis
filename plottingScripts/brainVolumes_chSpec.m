function [brainVolumes, mni_vox] = brainVolumes_chSpec(params, struct_plotBrain, channelStruct, fixedVols, inds_vol, groupInfo)
% updates/loads channel-specific brain volumes for slices plotting (defined for each channel separately)
% updates if fixedVols{vol}.loaded = 0
% input:
%   - params struct
%   - struct_plotBrain = params.plot_brainSlices OR params.plot_brain3D
%   - channelStruct = H.channels(ch) from header H
%   - fixedVols = pre-loaded brain volumes common for all channels (see: brainVolumes_fixed.m)
%   - inds_vol = volume indices to the common (catted) colormap (see: brainSlices_vols_fixed.m)
% output
%   - brainVolumes{vol} = 
%   struct with fields:
%             hdr: [1×1 struct]
%             vol: [177×209×176 double]
%             xyz: [3×6510768 double]
%              VI: [177×209×176 double]
%              xi: [1×177 double]
%              yi: [1×209 double]
%              zi: [1×176 double]
%     voxSize_new: 1
%          loaded: 1
%           cInds: [177×209×176 double]
%           aVals: [177×209×176 double]
%          vox_ix: 92
%          vox_iy: 120
%          vox_iz: 136
%   - mni_vox = MNI coordinates for the cross-hair position

% (c) Jiri, Nov21

%% set (or load) channel-specific brain volumes -> brainVolumes
nVolumes = size(struct_plotBrain.volumes2plot,1);
brainVolumes = cell(1, nVolumes);
for vol = 1:nVolumes
    if ~fixedVols{vol}.loaded
        % define the channel-specific brain volume (empty if not found)
%         file2load = getChSpecBrainVolName(channelStruct, struct_plotBrain.volumes2plot{vol,1}, groupInfo);  % TBD !!!
        brainAtlas_name = strrep(struct_plotBrain.volumes2plot{vol,1}, 'ch-spec_', '');
        file2load = ch2roi_assignAA(brainAtlas_name, channelStruct, groupInfo);

        % >>> load brain volume -> brain{vol} <<<
        if ~strcmp(file2load, 'n.a.')
            brainVolumes{vol} = brainVolumes_load(params, file2load, struct_plotBrain);
            V = brainVolumes{vol}.VI;
            
            % indices of color values for brain volume -> cVals
            brainVolumes{vol}.cInds = cVals2cInds(V, [prctile(V(:),1),prctile(V(:),99)], inds_vol(vol,:));

            % indices of transperency values for brain volume -> aVals
            brainVolumes{vol}.aVals = linTransform(V, [prctile(V(:),1),prctile(V(:),99.9)], struct_plotBrain.volumes2plot{vol,3});    
        end    
    else
        brainVolumes{vol} = fixedVols{vol};
    end
end

%% MNI coors (on the axis in brainVolumes{vol}.xi) -> mni_vox
mni_vox = [channelStruct.MNI_x, channelStruct.MNI_y, channelStruct.MNI_z];
for vol = 1:nVolumes
    if ~isempty(brainVolumes{vol}) 
        [vox_ix,vox_iy,vox_iz] = mni2vox(-mni_vox(1), mni_vox(2), mni_vox(3), brainVolumes{vol}.xi, brainVolumes{vol}.yi, brainVolumes{vol}.zi);
        brainVolumes{vol}.vox_ix = vox_ix;
        brainVolumes{vol}.vox_iy = vox_iy;
        brainVolumes{vol}.vox_iz = vox_iz;
    end
end
