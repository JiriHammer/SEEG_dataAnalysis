function [fixedVols, clrmap, inds_vol] = brainVolumes_fixed(params, struct_plotBrain)
% returns/loads pre-defined brain volumes (fixed, same for all channels) & colormap 
% with indices for each brain volume specied in:  struct_plotBrain.volumes2plot
% input:
%   - params
%   - struct_plotBrain = params.plot_brainSlices OR params.plot_brain3D
% fixedVols{vol} = 
%   struct with fields:
%             hdr: [1×1 struct]
%             vol: [253×298×251 double]
%             xyz: [3×18923894 double]
%              VI: [253×298×251 double]
%              xi: [1×253 double]
%              yi: [1×298 double]
%              zi: [1×251 double]
%     voxSize_new: 0.7000
%          loaded: 1
%           cInds: [253×298×251 double]
%           aVals: [253×298×251 double]
% clrmap = catted colormaps for each volume: N x 3 (RGB)
% inds_vol = nVolumes x 2 (indices to each volume colormap)
%
% Avoids thus loading the MRI or CT info for each channel again and again.
% If the brain volume is channel specific, returns: fixedVols{vol}.loaded = 0 

% (c) Jiri, Nov21

%% initialize
nVolumes = size(struct_plotBrain.volumes2plot,1);
params.plot_brainTopo.size_interpolate = 0.7;       % TO DO ...
clrmap = [];
inds_vol = [];
i_clr = 0;

%% load brain data -> fixedVols & clrmap
fixedVols = cell(1,nVolumes);    % pre-defined brain volumes (fixed, same for all channels)
for vol = 1:nVolumes
    
    % >>> load brain volume -> brain{vol} <<<
    file2load = struct_plotBrain.volumes2plot{vol,1};  % e.g.: wrCT, wT1, wT2, wc1T1, wc2T2, ...
    fixedVols{vol} = brainVolumes_load(params, file2load, struct_plotBrain);
    
    % colormap for brain volume -> clrmap
    this_clrmap = struct_plotBrain.volumes2plot{vol,2};
    clrmap = cat(1, clrmap,  this_clrmap);      % colormap for all brain volumes
    inds_vol = cat(1, inds_vol, [1, size(this_clrmap,1)]+i_clr);
    i_clr = i_clr + size(this_clrmap,1);
    
    % colormap for one brain volume
    fixedVols{vol}.clrmap = this_clrmap;        
    fixedVols{vol}.alfmap = struct_plotBrain.volumes2plot{vol,3};
    
    % map brain volumes to color (tranp.) indices based on their values (intensity)
    if fixedVols{vol}.loaded
        V = fixedVols{vol}.VI;
        
        % indices of color values for brain volume -> cVals
        fixedVols{vol}.cInds = cVals2cInds(V, [prctile(V(:),1),prctile(V(:),99)], inds_vol(end,:));

        % indices of transperency values for brain volume -> aVals
        fixedVols{vol}.aVals = linTransform(V, [prctile(V(:),1),prctile(V(:),99.9)], struct_plotBrain.volumes2plot{vol,3});    
    end
end
