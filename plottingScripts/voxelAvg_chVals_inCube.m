function [chVals_out, chMNIs_out] = voxelAvg_chVals_inCube(chVals, chnlsMNI, subj_inds, d_mm, n_subjContribute)
% finds average value from all channels that are located within a cube of d_mm

% (c) Jiri, Jan19

%% number of different subjects in each larger voxel
if nargin < 5
    n_subjContribute = 2;
end

%% define a larger voxel for smoothing (running averging)
if nargin < 4
    d_mm = 3;
end


%% sanity checks
assert(size(chnlsMNI,1) == size(chVals,1));
chnlsMNI = chnlsMNI';
assert(size(subj_inds,1) == size(chVals,1));

%% find "close-enough" channels in MNI & assign avg value
chVals_out = []; 
chMNIs_out = [];
for ch = 1:size(chVals,1)
    thisMNI = chnlsMNI(:,ch);
    i_chnls = find(abs(chnlsMNI(1,:) - thisMNI(1)) <= d_mm & ...
                   abs(chnlsMNI(2,:) - thisMNI(2)) <= d_mm & ...
                   abs(chnlsMNI(3,:) - thisMNI(3)) <= d_mm); 
    nSubj = numel(unique(subj_inds(i_chnls)));
    if nSubj >= n_subjContribute
        chVals_out = cat(1, chVals_out, mean(chVals(i_chnls),1));
        chMNIs_out = cat(2, chMNIs_out, chnlsMNI(:,ch));
    end
end

%% uncomment to debug
% figure;
% hold on;
% plot3(thisMNI(1), thisMNI(2), thisMNI(3), 'k*');
% for ch = 1:size(chVals,1)
%     plot3(chnlsMNI(1,ch), chnlsMNI(2,ch), chnlsMNI(3,ch), 'b.');
% end
% for ch = i_chnls
%     plot3(chnlsMNI(1,ch), chnlsMNI(2,ch), chnlsMNI(3,ch), 'ro');
% end
% axis equal

