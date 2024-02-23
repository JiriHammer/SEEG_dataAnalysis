function i_aa = assign_ch2aa(chnlStruct_H, cytoarchAreas, subjTag, chnlNum)
% converts channel topology to (larger) anatomical areas 'aa'
% input: channel info from header H = H.channel(N)
% output: vektor of AA indices

% (c) Jiri, Jun17

if nargin < 4, chnlNum = 0; end
if nargin < 3, subjTag = 'noName'; end

%% define AA = 'thisAA' (based on the latest Anatomy toolbox)
elAss = se_TabList_jiri(chnlStruct_H.MNI_x, chnlStruct_H.MNI_y, chnlStruct_H.MNI_z);
AA_H_struct = chnlStruct_H.ass_cytoarchMap;
AA_toolbox = elAss.probabAss_area{1}; 
if ~strcmp(AA_H_struct, elAss.probabAss_area{1})
    if strcmp(AA_H_struct, 'aufCS') || strcmp(AA_H_struct, 'auf CS')
        thisAA = 'Motor';
    elseif strcmp(AA_toolbox, 'n.a.')
        thisAA = getLargerAnatomicalArea(AA_H_struct);
    else
        thisAA = getLargerAnatomicalArea(AA_toolbox);
    end
else
    thisAA = getLargerAnatomicalArea(AA_toolbox);
end

%% check if it was assigned?
[aa_found,i_aa] = ismember(thisAA, cytoarchAreas);
if ~aa_found
    display([subjTag ', ch = ' num2str(chnlNum) ', not assigned: H struct = ' chnlStruct_H.ass_cytoarchMap ', toolbox = ' elAss.probabAss_area{1}]);
    i_aa = nan;
end
