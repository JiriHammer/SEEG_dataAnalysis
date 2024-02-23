function [chVals_out, chMNIs_out] = voxelAvg_chVals_3Dconv(chVals, plotInfo, d_mm)
% returns a running average of channels values in MNI space
% define a larger voxel
% assign voxels by channels values (chVals)
% smooth/average (3D convolution) & division by the number of times the values
% were assigned in each voxel

% (c) Jiri, Nov18

%% define a larger voxel for smoothing (running averging)
if nargin < 3
    d_mm = 3;
end

%% defaults
if ~isfield(plotInfo, 'subj_inds') 
    plotInfo.subj_inds = ones(1,size(plotInfo.chnlsMNI,2));
    n_subjContribute = 1;
else
    n_subjContribute = 2;
end

%% sanity checks
assert(size(plotInfo.chnlsMNI,2) == size(chVals,1));
assert(size(plotInfo.subj_inds,2) == size(chVals,1));

%% assign voxels in MNI space
V_vals = zeros(size(plotInfo.brain.VI));    % values for running avg (convolution)
V_count = zeros(size(plotInfo.brain.VI));    % values counter (if same voxel has two values from diff. subjects)
V_used = zeros(size(plotInfo.brain.VI));    % values used
V_subj = zeros(size(plotInfo.brain.VI));    % contributions from different subjects
for ch = 1:size(chVals,1)
    ix = closestval(plotInfo.brain.xi, plotInfo.chnlsMNI(1,ch));
    iy = closestval(plotInfo.brain.yi, plotInfo.chnlsMNI(2,ch));
    iz = closestval(plotInfo.brain.zi, plotInfo.chnlsMNI(3,ch));
    V_used(ix,iy,iz) = 1;
    V_count(ix,iy,iz) = V_count(ix,iy,iz)+1;
    V_vals(ix,iy,iz) = V_vals(ix,iy,iz)+chVals(ch);
end

%% convolve: running sum in a larger (smoothing) voxel
d_vox = closestval(plotInfo.brain.xi,d_mm) - closestval(plotInfo.brain.xi,0.0); % assumes same axis steps
C_conv = ones(d_vox,d_vox,d_vox);   % convolution kernel: larger voxel for convolution
V_vals_conv = convn(V_vals,C_conv,'same');
V_used_conv = convn(V_used,C_conv,'same');
V_count_conv =convn(V_count,C_conv,'same');

%% avg over used values
I = V_used_conv ~= 0;                           % used indices (used at least once)
%I = V_used_conv >= 2;                           % used indices (used in more than 2 channels)
V_conv = zeros(size(plotInfo.brain.VI));        % init
V_conv(I) = V_vals_conv(I)./V_count_conv(I);    % avg

%% individual subject contributions
V_subj = zeros(size(plotInfo.brain.VI));    % contributions from different subjects
i_subj_unq = unique(plotInfo.subj_inds);
for subj = 1:size(i_subj_unq,2)
    thisSubj = i_subj_unq(subj);
    ch_thisSubj = find(plotInfo.subj_inds == thisSubj);
    
    % assign 1 to MNI coors ---> V_tmp
    V_tmp = zeros(size(plotInfo.brain.VI));    % contributions from different subjects
    for ch = ch_thisSubj
        ix = closestval(plotInfo.brain.xi, plotInfo.chnlsMNI(1,ch));
        iy = closestval(plotInfo.brain.yi, plotInfo.chnlsMNI(2,ch));
        iz = closestval(plotInfo.brain.zi, plotInfo.chnlsMNI(3,ch));
        V_tmp(ix,iy,iz) = 1;
    end
    
    % convolve (smooth over volume) ---> V_tmp_conv
    V_tmp_conv = convn(V_tmp,C_conv,'same');
    
    % used indices (used at least once) in larger voxel
    V_tmp_conv(V_tmp_conv ~= 0) = 1;
    
    % add to all subjects used
    V_subj = V_subj + V_tmp;
end

%% return channel values at MNI voxels
chVals_out = [];
chMNIs_out = [];
for ch = 1:size(chVals,1)
    ix = closestval(plotInfo.brain.xi, plotInfo.chnlsMNI(1,ch));
    iy = closestval(plotInfo.brain.yi, plotInfo.chnlsMNI(2,ch));
    iz = closestval(plotInfo.brain.zi, plotInfo.chnlsMNI(3,ch));
    assert(V_used_conv(ix,iy,iz) > 0);
    assert(V_count_conv(ix,iy,iz) > 0);
    if V_subj(ix,iy,iz) >= n_subjContribute
        chVals_out = cat(1, chVals_out, V_conv(ix,iy,iz));
        chMNIs_out = cat(2, chMNIs_out, plotInfo.chnlsMNI(:,ch));
    else
        1+1;
    end
end

%% debug figure
% z_slice = 65;
% figure;
% subplot(1,2,1); hold on;
% imagesc(plotInfo.brain.yi, plotInfo.brain.xi, V_vals(:,:,z_slice));
% axis tight;
% axis equal;
% 
% subplot(1,2,2); hold on;
% imagesc(plotInfo.brain.yi, plotInfo.brain.xi, V_conv(:,:,z_slice));
% axis tight;
% axis equal;

