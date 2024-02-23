function AA_list = anatomicalAreas_getList(listName, groupInfo)
% returns list of anatomical areas based on specification in 'listName'
% each row in the cell array 'AA_list' is 1 anatomical area (AA)
% first column of 'AA_list' is the name of the area
% second column of 'AA_list' is a cell array of strings defining the area
% for example:
%    'Insula'            {1x3  cell}, where
%       {1x3  cell} = 'Area Id1 '    'Area Ig2 '    'Area Ig1 '

% (c) Jiri, Jan18

if nargin < 2
    groupInfo = struct;
end

%% anatomy assignments from neurologists (Adam Kalina) -> must contain 'LabelsAK'
if ~isempty(strfind(listName, 'LabelsAK'))
    fileName = groupInfo.file_brainClusters;
    assert(exist(fileName,'file') == 2);  
    [numdata, txtdata, alldata] = xlsread(fileName);
    % remove spaces at the end of each line
    for aa = 1:size(alldata,1)
        if strcmp(alldata{aa}(end),' ')
            alldata{aa} = alldata{aa}(1:end-1);
        end
    end
    
    % make list of AAs
    AA_list = cell(size(alldata,1),2);
    for aa = 1:size(alldata,1)
        AA_list{aa,1} = alldata{aa,1};
        for n = 2:size(alldata,2)
            if ~isnan(alldata{aa,n})
                AA_list{aa,2}{1,n-1} = str_splitSortDelWM(alldata{aa,n});   % sort & remove wm
            end
        end
        AA_list{aa,2} = unique(lower(AA_list{aa,2}));
    end
end

%% Yeo7 networks
% if ~isempty(strfind(listName, 'Yeo7')) % strcmp(listName, 'Yeo7')
%     AA_list = {
%     'Visual', {'Visual'};
%     'Somatomotor', {'Somatomotor'};
%     'Dorsal Attention', {'Dorsal Attention'};
%     'Ventral Attention', {'Ventral Attention'};
%     'Limbic', {'Limbic'};
%     'Frontoparietal', {'Frontoparietal'};
%     'Default', {'Default'};    
%     };
% end

%% TBD: Yeo7 networks (re-defined order for MK's EAN poster)
% if ~isempty(strfind(listName, 'Yeo7')) % strcmp(listName, 'Yeo7')
%     AA_list = {
%     'Default', {'Default'};      
%     'Dorsal Attention', {'Dorsal Attention'};
%     'Ventral Attention', {'Ventral Attention'};
%     'Frontoparietal', {'Frontoparietal'};
%     'Somatomotor', {'Somatomotor'};
%     'Limbic', {'Limbic'};
%     'Visual', {'Visual'};
%     };
% end

%% Yeo7 networks (order based on spectra similarity in SEI ms)
if ~isempty(strfind(listName, 'Yeo7')) % strcmp(listName, 'Yeo7')
    AA_list = {
    'Visual', {'Visual'};
    'Dorsal Attention', {'Dorsal Attention'};
    'Somatomotor', {'Somatomotor'};    
    'Ventral Attention', {'Ventral Attention'};
    'Frontoparietal', {'Frontoparietal'};
    'Limbic', {'Limbic'};
    'Default', {'Default'};    
    };
end

%% Yeo17 networks
% if ~isempty(strfind(listName, 'Yeo17')) % strcmp(listName, 'Yeo17')
%     AA_list = {
%     'Visual peripheral', {'Visual peripheral'};
%     'Visual central', {'Visual central'};
%     'Somatomotor A', {'Somatomotor A'};
%     'Somatomotor B', {'Somatomotor B'};
%     'Dorsal attention A', {'Dorsal attention A'};
%     'Dorsal attention B', {'Dorsal attention B'};
%     'Ventral attention', {'Ventral attention'};
%     'Salience', {'Salience'};
%     'Limbic N9', {'Limbic N9'};
%     'Limbic N10', {'Limbic N10'};
%     'Control C', {'Control C'};
%     'Control A', {'Control A'};
%     'Control B', {'Control B'};
%     'Default D (Auditory)', {'Default D (Auditory)'};
%     'Default C', {'Default C'};
%     'Default A', {'Default A'};
%     'Default B', {'Default B'};
%     };
% end

%% Yeo17 networks (order based on spectra similarity in SEI ms)
if ~isempty(strfind(listName, 'Yeo17')) % strcmp(listName, 'Yeo17')
    AA_list = {
    'Visual peripheral', {'Visual peripheral'};
    'Visual central', {'Visual central'};
    'Dorsal attention A', {'Dorsal attention A'};
    'Dorsal attention B', {'Dorsal attention B'};    
    'Somatomotor A', {'Somatomotor A'};
    'Somatomotor B', {'Somatomotor B'};
    'Ventral attention', {'Ventral attention'};
    'Salience', {'Salience'};
    'Control C', {'Control C'};
    'Control A', {'Control A'};
    'Control B', {'Control B'};
    'Limbic N9', {'Limbic N9'};
    'Limbic N10', {'Limbic N10'};    
    'Default D (Auditory)', {'Default D (Auditory)'};
    'Default C', {'Default C'};
    'Default A', {'Default A'};
    'Default B', {'Default B'};
    };
end

%% Destrieux atlas (wo laterality)
if ~isempty(strfind(listName, 'Destrie')) % strcmp(listName, 'Destrie')
    if ~isfield(groupInfo, 'path_isarg_atlas')
        params_default;
        groupInfo.path_isarg_atlas = params.path2others.isarg_atlas;
    end
    file_roiTable  = [groupInfo.path_isarg_atlas filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep 'destrieux2009_rois_labels.csv'];    
    assert(exist(file_roiTable,'file') == 2);    
    T_roi = readtable(file_roiTable);      
    
    AA_list = cell(size(T_roi,1),2);
    for aa = 1:size(T_roi,1)
        AA_list{aa,1} = T_roi.name{aa};     % string
        AA_list{aa,2} = T_roi.name(aa);     % cell
    end
end

%% Destrieux atlas (wi laterality)
if ~isempty(strfind(listName, 'DestLat')) % strcmp(listName, 'DestLat')
    if ~isfield(groupInfo, 'path_isarg_atlas')
        params_default;
        groupInfo.path_isarg_atlas = params.path2others.isarg_atlas;
    end    
    file_roiTable  = [groupInfo.path_isarg_atlas filesep 'coregistration_isarg\atlas_tmp\Destrieux_mni152' filesep 'destrieux2009_rois_labels_lateralized.csv'];    
    assert(exist(file_roiTable,'file') == 2);    
    T_roi = readtable(file_roiTable);      
    
    AA_list = cell(size(T_roi,1),2);
    for aa = 1:size(T_roi,1)
        AA_list{aa,1} = T_roi.name{aa};     % string
        AA_list{aa,2} = T_roi.name(aa);     % cell
    end
end

%% Mars atlas (wo laterality) (note: for lateralized version, use Mars_addLat)
if ~isempty(strfind(listName, 'Mars')) 
    if ~isfield(groupInfo, 'path_isarg_atlas')
        params_default;
        groupInfo.path_isarg_atlas = params.path2others.isarg_atlas;
    end
    file_roiTable  = [groupInfo.path_isarg_atlas filesep 'coregistration_isarg\atlas_tmp\Mars_mni152' filesep 'Mars_labels.csv'];    
    assert(exist(file_roiTable,'file') == 2);    
    T_roi = readtable(file_roiTable);      
    
    AA_list = cell(size(T_roi,1),2);
    for aa = 1:size(T_roi,1)
        AA_list{aa,1} = T_roi.Label{aa};     % string
        AA_list{aa,2} = T_roi.Label(aa);     % cell
    end
end

%% MARS-msDist = user modified Mars atlas: for msDist
if ~isempty(strfind(listName, 'MARS-msDist')) 
    AA_list = {
    'Premotor', {'PMdm','PMdl','PFcdm'};
    'Motor', {'Mdl','Mdm'};
    'Somatosensory', {'Sdl','Sdm'};
    'IPC', {'IPCd','IPCv'};
    'SPC', {'SPC','SPCm'};
...    'Insula', {'Insula'};
    };
end
    
%% add laterality (L/R) ?
assert(exist('AA_list','var') == 1);
if ~isempty(strfind(listName, '_addLat'))
    AA_tmp = AA_list;
    % make list of AAs
    AA_list = cell(2*size(AA_tmp,1),2);
    n = 1;
    for aa = 1:size(AA_tmp,1)
        % Left hemisphere
        AA_list{n,1} = ['L ' AA_tmp{aa,1}];
        for c = 1:size(AA_tmp{aa,2},2)
            AA_list{n,2}{1,c} = ['L ' AA_tmp{aa,2}{1,c}];
        end
        n = n+1;
        
        % Right hemisphere
        AA_list{n,1} = ['R ' AA_tmp{aa,1}];
        for c = 1:size(AA_tmp{aa,2},2)
            AA_list{n,2}{1,c} = ['R ' AA_tmp{aa,2}{1,c}];
        end
        n = n+1;
    end
%     for aa = 1:size(AA_tmp,1)
%         AA_list{n,1} = ['R ' AA_tmp{aa,1}];
%         for c = 1:size(AA_tmp{aa,2},2)
%             AA_list{n,2}{1,c} = ['R ' AA_tmp{aa,2}{1,c}];
%         end
%         n = n+1;
%     end
end

return;











%% OLD & TO DO (now = 01.08.2022)
%% cyto-architectonic assignments (motor plus sensory area combined)
if ~isempty(strfind(listName, 'cytoarch_large_motorPlusSensory')) % strcmp(listName, 'cytoarch_large_motorPlusSensory')
    AA_list = {
    'FrontalPole',  {'Area Fp1','Area Fp2'};
    'BasalForebrain',{'BF (Ch 1-3)','BF (Ch 4)'};
    'OFC',          {'Area Fo1', 'Area Fo2','Area Fo3'};
    'Broca',        {'Area 44','Area 45'};
    'MotorSensory', {'Area 4a','Area 4p',  ...                  % M1 (v22c)
                     'Area 6','aufCS','auf CS',...              % M1 (v18 & FR hierarchical assignment)
                     'Area 1','Area 2','Area 3a','Area 3b'};    % S1 
    'Insula',       {'Area Id1 ','Area Ig2 ','Area Ig1 '};
    'Cingulum',     {'Area 25','Area s24','Area s32','Area 33'};
    'Operculum',    {'Area OP1 [SII]','Area OP2 [PIVC]','Area OP3 [VS]','Area OP4 [PV]'};
    'Auditory',     {'Area TE 1.0','Area TE 1.1','Area TE 1.2','Area TE 3',    'TE 3'};
    'Amygdala',     {'Amygdala (CM)','Amygdala (SF)','Amygdala (AStr)','Amygdala (LB)'};
    'Hippocampus',  {'DG (Hippocampus)','CA1 (Hippocampus)','CA2 (Hippocampus)','CA3 (Hippocampus)','HATA Region','Entorhinal Cortex','Subiculum'};
    'IPL',          {'Area PFop (IPL)','Area PFm (IPL)','Area PF (IPL)','Area PFcm (IPL)','Area PGp (IPL)','Area PGa (IPL)','Area PFt (IPL)',     'IPC (PF)'};
    'SPL/IPS',      {'Area 5M (SPL)','Area 5L (SPL)','Area 5Ci (SPL)','Area 7M (SPL)','Area 7A (SPL)','Area 7PC (SPL)','Area 7P (SPL)','SPL (7PC)',   'Area hIP1 (IPS)','Area hIP3 (IPS)','Area hIP2 (IPS)'};
    'Visual',       {'Area FG1','Area FG2','Area  FG3','Area FG4','Area hOc1 [V1]','Area hOc2 [V2]','Area hOc4v [V4(v)]', 'Area hOc5 [V5/MT]','Area hOc3v [V3v]','Area hOc4d [V3A]','Area hOc4la','Area hOc3d [V3d]','Area hOc4lp'};
    'Thalamus',     {'Thal: Prefrontal','Thal: Motor','Thal: Parietal','Thal: Premotor','Thal: Somatosensory','Thal: Visual','Thal: Temporal'};
    'Cerebellum',   {'Lobule IX (Hem)','Lobule VI (Verm)','Lobule VIIa crusI (Verm)','Lobule VIIa crusII (Verm)','Lobule VIIIb (Verm)','Lobule VIIb (Hem)','Lobule VIIIa (Verm)','Lobule X (Verm)','Lobule VIIIb (Hem)','Lobule VIIb (Verm)','Lobule VIIIa (Hem)','Lobule IX (Verm)','Lobule VIIa crusI (Hem)','Lobule I IV (Hem)','Lobule VIIa crusII (Hem)','Lobule X (Hem)','Lobule VI (Hem)','Lobule V (Hem)','Ventral Dentate Nucleus','Dorsal Dentate Nucleus','Interposed Nucleus','Fastigii Nucleus'};
    };
end

%% cyto-architectonic assignments
if ~isempty(strfind(listName, 'cytoarch_large')) % strcmp(listName, 'cytoarch_large')
    AA_list = {
    'FrontalPole',  {'Area Fp1','Area Fp2'};
    'BasalForebrain',{'BF (Ch 1-3)','BF (Ch 4)'};
    'OFC',          {'Area Fo1', 'Area Fo2','Area Fo3'};
    'Broca',        {'Area 44','Area 45'};
    'Motor',        {'Area 4a','Area 4p',    'Area 6'};
    'Sensory',      {'Area 1','Area 2','Area 3a','Area 3b'};
    'Insula',       {'Area Id1 ','Area Ig2 ','Area Ig1 '};
    'Cingulum',     {'Area 25','Area s24','Area s32','Area 33'};
    'Operculum',    {'Area OP1 [SII]','Area OP2 [PIVC]','Area OP3 [VS]','Area OP4 [PV]'};
    'Auditory',     {'Area TE 1.0','Area TE 1.1','Area TE 1.2','Area TE 3',    'TE 3'};
    'Amygdala',     {'Amygdala (CM)','Amygdala (SF)','Amygdala (AStr)','Amygdala (LB)'};
    'Hippocampus',  {'DG (Hippocampus)','CA1 (Hippocampus)','CA2 (Hippocampus)','CA3 (Hippocampus)','HATA Region','Entorhinal Cortex','Subiculum'};
    'IPL',          {'Area PFop (IPL)','Area PFm (IPL)','Area PF (IPL)','Area PFcm (IPL)','Area PGp (IPL)','Area PGa (IPL)','Area PFt (IPL)',     'IPC (PF)'};
    'SPL/IPS',      {'Area 5M (SPL)','Area 5L (SPL)','Area 5Ci (SPL)','Area 7M (SPL)','Area 7A (SPL)','Area 7PC (SPL)','Area 7P (SPL)','SPL (7PC)',   'Area hIP1 (IPS)','Area hIP3 (IPS)','Area hIP2 (IPS)'};
    'Visual',       {'Area FG1','Area FG2','Area  FG3','Area FG4','Area hOc1 [V1]','Area hOc2 [V2]','Area hOc4v [V4(v)]', 'Area hOc5 [V5/MT]','Area hOc3v [V3v]','Area hOc4d [V3A]','Area hOc4la','Area hOc3d [V3d]','Area hOc4lp'};
    'Thalamus',     {'Thal: Prefrontal','Thal: Motor','Thal: Parietal','Thal: Premotor','Thal: Somatosensory','Thal: Visual','Thal: Temporal'};
    'Cerebellum',   {'Lobule IX (Hem)','Lobule VI (Verm)','Lobule VIIa crusI (Verm)','Lobule VIIa crusII (Verm)','Lobule VIIIb (Verm)','Lobule VIIb (Hem)','Lobule VIIIa (Verm)','Lobule X (Verm)','Lobule VIIIb (Hem)','Lobule VIIb (Verm)','Lobule VIIIa (Hem)','Lobule IX (Verm)','Lobule VIIa crusI (Hem)','Lobule I IV (Hem)','Lobule VIIa crusII (Hem)','Lobule X (Hem)','Lobule VI (Hem)','Lobule V (Hem)','Ventral Dentate Nucleus','Dorsal Dentate Nucleus','Interposed Nucleus','Fastigii Nucleus'};
    };
end

%% TO DO: cyto-architectonic assignments of individual areas in Anatomy toolbox
