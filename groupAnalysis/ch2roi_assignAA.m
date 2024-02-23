function thisAA = ch2roi_assignAA(brainAtlas_name, chAss_isargAtlas, groupInfo, list_AA)
% returns thisAA (string) of assigned anatomical area (AA aka ROI) 
% thisAA = a row in the first column of 'list_AA' or 'n.a.' (if not assigned)
% based on 'brainAtlas_name' (case insensitive!), selects the corresponding field in
% 'chAss_isargAtlas' struct and looks for case insensitive string comparison (ROI + distance tolerance)
% for details, see: anatomicalArea_getName.m
% also adds lateralization (prefix = 'L ' or 'R '), if the 'brainAtlas_name' contains '_addLat'

% (c) Jiri, Aug22

if nargin < 4
    list_AA = anatomicalAreas_getList(brainAtlas_name, groupInfo);   % define ROIs (anatom. areas) 
end

%% max distance tolerance to ROI
if ~isfield(groupInfo, 'maxDistToROI')
    groupInfo.maxDistToROI = 0;
end

%% add lateralization
addLateralization = false;
if ~isempty(strfind(brainAtlas_name, '_addLat'))
    addLateralization = true;
end

%% find the right atlas and assign -> thisAA
if ~isempty(strfind(lower(brainAtlas_name), 'yeo7')) % strcmp(brainAtlas_name,'Yeo7')   
    thisAA = anatomicalArea_getName('yeo7_name', list_AA, chAss_isargAtlas, addLateralization, groupInfo.maxDistToROI);
    
elseif ~isempty(strfind(lower(brainAtlas_name), 'yeo17')) % strcmp(brainAtlas_name,'Yeo17')   
    thisAA = anatomicalArea_getName('yeo17_name', list_AA, chAss_isargAtlas, addLateralization, groupInfo.maxDistToROI);
    
elseif ~isempty(strfind(lower(brainAtlas_name), 'destrie')) % strcmp(brainAtlas_name,'Destrie')   
    thisAA = anatomicalArea_getName('destrie_name', list_AA, chAss_isargAtlas, addLateralization, groupInfo.maxDistToROI);
    
elseif ~isempty(strfind(lower(brainAtlas_name), 'destlat')) % strcmp(brainAtlas_name,'DestLat')   
    thisAA = anatomicalArea_getName('destLat_name', list_AA, chAss_isargAtlas, addLateralization, groupInfo.maxDistToROI);
    
elseif ~isempty(strfind(lower(brainAtlas_name), 'mars')) % note: for lateralized version, use Mars_addLat
    thisAA = anatomicalArea_getName('mars_name', list_AA, chAss_isargAtlas, addLateralization, groupInfo.maxDistToROI);
    
elseif ~isempty(strfind(lower(brainAtlas_name), 'labelsak')) % strcmp(brainAtlas_name,'LabelsAK')   
    thisAA = anatomicalArea_getName('LabelsAK', list_AA, chAss_isargAtlas, addLateralization, groupInfo.maxDistToROI);
    
else
    disp(['WARNING: atlas = ' brainAtlas_name ' not found!']);
    thisAA = 'n.a.';
end
