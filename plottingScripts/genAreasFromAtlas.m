%% generating script for the areas in each atlas
%   -> Yeo7
%   -> Yeo17
%   -> Destrieux (incl. laterality)
%   -> Mars (incl. laterality)

% (c) Jiri, Aug22

%% Yeo7 Networks
outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\yeo7';
assert(exist(outDir,'dir') == 7);
cd(outDir);

% load .nii file
pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v7';
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
            'wYeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
assert(exist(fileName,'file') == 2);
V=spm_vol(fileName);
VOL=spm_read_vols(V);

% load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
            'wYeo2011_7Networks.mat'];
assert(exist(fileName,'file') == 2);
load(fileName,'Yeo_7N_labels');
AA_inds = Yeo_7N_labels;

% separate into NII files for each AA
for i=1:length(AA_inds)
    VN=V;
    VN.fname=[AA_inds{i} '.nii'];
    
    VOLN=VOL==i;
    
    spm_create_vol(VN);
    spm_write_vol(VN,VOLN);
end
disp(['Yeo7 areas (.nii volumes) saved in: ' outDir]);

%% Yeo17 Networks
pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v7';
outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\yeo17';
assert(exist(outDir,'dir') == 7);
cd(outDir);

% load .nii file
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
            'wYeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
assert(exist(fileName,'file') == 2);
V=spm_vol(fileName);
VOL=spm_read_vols(V);

% load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
            'wYeo2011_17Networks.mat'];
assert(exist(fileName,'file') == 2);
load(fileName,'Yeo_17N_labels');
AA_inds = Yeo_17N_labels;

% separate into NII files for each AA
for i=1:length(AA_inds)
    VN=V;
    VN.fname=[AA_inds{i} '.nii'];
    
    VOLN=VOL==i;
    
    spm_create_vol(VN);
    spm_write_vol(VN,VOLN);
end
disp(['Yeo17 areas (.nii volumes) saved in: ' outDir]);

%% Destrieux atlas
pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v7';
outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\destrie';
assert(exist(outDir,'dir') == 7);
cd(outDir);

% load .nii file
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep ...
            'wdestrieux2009_rois.nii'];
assert(exist(fileName,'file') == 2);
V=spm_vol(fileName);
VOL=spm_read_vols(V);

% load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
file_roiTable  = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep ...
                 'destrieux2009_rois_labels.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable); 
AA_name = cell(size(T_roi,1),1);
AA_inds = nan(size(T_roi,1),1);
for aa = 1:size(T_roi,1)
    AA_name{aa,1} = T_roi.name{aa};
    AA_inds(aa,1) = T_roi.index(aa); 
end

% separate into NII files for each AA
for i=2:size(AA_inds,1)     % starts with 2 (1 = Unknown)
    VN=V;
    VN.fname=[AA_name{i} '.nii'];
    
    VOLN=VOL==AA_inds(i);
    
    spm_create_vol(VN);
    spm_write_vol(VN,VOLN);
end
disp(['Destrieux areas (.nii volumes) saved in: ' outDir]);

%% Destrieux lateralized atlas
pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v7';
outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\destLat';
assert(exist(outDir,'dir') == 7);
cd(outDir);

% load .nii file
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep ...
            'wdestrieux2009_roisi_lateralized.nii'];
assert(exist(fileName,'file') == 2);
V=spm_vol(fileName);
VOL=spm_read_vols(V);

% load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
file_roiTable  = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep ...
                 'destrieux2009_rois_labels_lateralized.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable); 
AA_name = cell(size(T_roi,1),1);
AA_inds = nan(size(T_roi,1),1);
for aa = 1:size(T_roi,1)
    AA_name{aa,1} = T_roi.name{aa};
    AA_inds(aa,1) = T_roi.index(aa); 
end

% separate into NII files for each AA
for i=2:size(AA_inds,1)     % starts with 2 (1 = Unknown)
    VN=V;
    VN.fname=[AA_name{i} '.nii'];
    
    VOLN=VOL==AA_inds(i);
    
    spm_create_vol(VN);
    spm_write_vol(VN,VOLN);
end
disp(['Destrieux lateralized areas (.nii volumes) saved in: ' outDir]);

%% Mars lateralized atlas
pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v8';
outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\marsLat';
assert(exist(outDir,'dir') == 7);
cd(outDir);

% load .nii file
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Mars_mni152' filesep ...
            'wcolin27_MNI_MarsAtlas_v2.nii'];
assert(exist(fileName,'file') == 2);
V=spm_vol(fileName);
VOL=spm_read_vols(V);

% load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
file_roiTable  = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Mars_mni152' filesep ...
                 'Mars_labels.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable); 
AA_name = [];
AA_inds = [];
n = 1;
for aa = 1:size(T_roi.LeftIndex,1)          % L hemisphere
    AA_name{n,1} = ['L ' T_roi.Label{aa}];
    AA_inds(n,1) = T_roi.LeftIndex(aa);
    n = n+1;
end
for aa = 1:size(T_roi.RightIndex,1)          % R hemisphere
    AA_name{n,1} = ['R ' T_roi.Label{aa}];
    AA_inds(n,1) = T_roi.RightIndex(aa);
    n = n+1;
end

% separate into NII files for each AA
for i=1:size(AA_inds,1)     % starts with 2 (1 = Unknown)
    VN=V;
    VN.fname=[AA_name{i} '.nii'];
    
    VOLN=VOL==AA_inds(i);
    
    spm_create_vol(VN);
    spm_write_vol(VN,VOLN);
end
disp(['Mars lateralized areas (.nii volumes) saved in: ' outDir]);

%% Mars atlas (no laterality)
pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v8';
outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\mars';
assert(exist(outDir,'dir') == 7);
cd(outDir);

% load .nii file
fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Mars_mni152' filesep ...
            'wcolin27_MNI_MarsAtlas_v2.nii'];
assert(exist(fileName,'file') == 2);
V=spm_vol(fileName);
VOL=spm_read_vols(V);

% load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
file_roiTable  = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Mars_mni152' filesep ...
                 'Mars_labels.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable); 
AA_name = [];
AA_inds = [];
n = 1;
for aa = 1:size(T_roi.LeftIndex,1)
    AA_name{n,1} = T_roi.Label{aa};
    AA_inds = cat(1, AA_inds, [T_roi.LeftIndex(aa), T_roi.RightIndex(aa)]);
    n = n+1;
end

% separate into NII files for each AA
for i=1:size(AA_inds,1)     % starts with 2 (1 = Unknown)
    VN=V;
    VN.fname=[AA_name{i} '.nii'];
    
    VOLN= (VOL==AA_inds(i,1)) | (VOL==AA_inds(i,2));    % join both areas
    
    spm_create_vol(VN);
    spm_write_vol(VN,VOLN);
end
disp(['Mars areas (.nii volumes) saved in: ' outDir]);








%% OLD (same as above) - del?
% %% generating script for the areas in each atlas
% %   -> Yeo7
% %   -> Yeo17
% 
% % (c) Jiri, Aug22
% 
% %% Yeo7 Networks
% outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\yeo7';
% assert(exist(outDir,'dir') == 7);
% cd(outDir);
% 
% % load .nii file
% pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v7';
% fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
%             'wYeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
% assert(exist(fileName,'file') == 2);
% V=spm_vol(fileName);
% VOL=spm_read_vols(V);
% 
% % load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
% fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
%             'wYeo2011_7Networks.mat'];
% assert(exist(fileName,'file') == 2);
% load(fileName,'Yeo_7N_labels');
% AA_inds = Yeo_7N_labels;
% 
% % separate into NII files for each AA
% for i=1:length(AA_inds)
%     VN=V;
%     VN.fname=[AA_inds{i} '.nii'];
%     
%     VOLN=VOL==i;
%     
%     spm_create_vol(VN);
%     spm_write_vol(VN,VOLN);
% end
% disp(['Yeo7 areas .nii volumes saved in: ' outDir]);
% 
% %% Yeo17 Networks
% pathBeg = 'G:\shareData\visualization_normBrains\isarg_atlas_v7';
% outDir = 'G:\shareData\visualization_normBrains\atlasAreas_MNI\yeo17';
% assert(exist(outDir,'dir') == 7);
% cd(outDir);
% 
% % load .nii file
% fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
%             'wYeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'];
% assert(exist(fileName,'file') == 2);
% V=spm_vol(fileName);
% VOL=spm_read_vols(V);
% 
% % load corresponding labels (!!! labels MUST correspond to indices in NII file !!!)
% fileName = [pathBeg filesep 'coregistration_isarg\atlas_tmp\Yeo_fmri_mni152' filesep ...
%             'wYeo2011_17Networks.mat'];
% assert(exist(fileName,'file') == 2);
% load(fileName,'Yeo_17N_labels');
% AA_inds = Yeo_17N_labels;
% 
% % separate into NII files for each AA
% for i=1:length(AA_inds)
%     VN=V;
%     VN.fname=[AA_inds{i} '.nii'];
%     
%     VOLN=VOL==i;
%     
%     spm_create_vol(VN);
%     spm_write_vol(VN,VOLN);
% end
% disp(['Yeo17 areas .nii volumes saved in: ' outDir]);
