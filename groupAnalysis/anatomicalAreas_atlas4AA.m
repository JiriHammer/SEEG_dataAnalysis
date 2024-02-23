function atlasName = anatomicalAreas_atlas4AA(aaName)
% returns name of the atlas given an anatomical area (or n.a. if not found)

% (c) Jiri, Aug22

atlasName = 'n.a.';    % init (not found)

%% was aaName part of the Yeo7 atlas?
AA_list = {
    'Visual';
    'Somatomotor';
    'Dorsal Attention';
    'Ventral Attention';
    'Limbic';
    'Frontoparietal';
    'Default'
};
if ismember(aaName, AA_list)
    atlasName = 'yeo7';
    return;
end

%% was aaName part of the Yeo17 atlas?
AA_list = {
    'Visual peripheral';
    'Visual central';
    'Somatomotor A';
    'Somatomotor B';
    'Dorsal attention A';
    'Dorsal attention B';
    'Ventral attention';
    'Salience';
    'Limbic N9';
    'Limbic N10';
    'Control C';
    'Control A';
    'Control B';
    'Default D (Auditory)';
    'Default C';
    'Default A';
    'Default B';
};
if ismember(aaName, AA_list)
    atlasName = 'yeo17';
    return;
end

%%  was aaName part of the Destrieux atlas (wo laterality)
params_default;
file_roiTable  = [params.path2others.isarg_atlas filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep 'destrieux2009_rois_labels.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable);      
AA_list = cell(size(T_roi,1),1);
for aa = 1:size(T_roi,1)
    AA_list{aa,1} = T_roi.name{aa};     % string
end

if ismember(aaName, AA_list)
    atlasName = 'destrie';
    return;
end

%% Destrieux atlas (wi laterality)
params_default;  
file_roiTable  = [params.path2others.isarg_atlas filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep 'destrieux2009_rois_labels_lateralized.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable);      
AA_list = cell(size(T_roi,1),1);
for aa = 1:size(T_roi,1)
    AA_list{aa,1} = T_roi.name{aa};     % string
end

if ismember(aaName, AA_list)
    atlasName = 'destLat';
    return;
end

%% Mars atlas (wo laterality) (note: for lateralized version, use Mars_addLat)
params_default;
file_roiTable  = [params.path2others.isarg_atlas filesep 'coregistration_isarg\atlas_tmp\Mars_mni152' filesep 'Mars_labels.csv'];    
assert(exist(file_roiTable,'file') == 2);    
T_roi = readtable(file_roiTable);      

AA_list = cell(size(T_roi,1),1);
for aa = 1:size(T_roi,1)
    AA_list{aa,1} = T_roi.Label{aa};     % string
end

if ismember(aaName, AA_list)
    atlasName = 'mars';
    return;
end

