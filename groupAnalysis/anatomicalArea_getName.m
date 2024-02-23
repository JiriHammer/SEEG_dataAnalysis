function aaName = anatomicalArea_getName(atlasName, AA_list, H_ch, addLateralization, maxDist)
% based on atlasName & assignments in header H, returns name of the assigned anatomical area aaName
% works for multiple atlases (Yeo, Destrieux, neurologists labels = LabelsAK)
% assignment is NOT case sensitive (Visual = visual, STG = stg)
% implemetation:
% (1) atlasName -> chAss (channels assignment in H.channels)
% (2) in case of "multiple area" assignments (typically from LabelsAK): 
%   - separate (separator = '/') to cell arrays 
%   - sort alphabetically
%   - remove 'wm'
% (3) - adds laterality (if selected)
% (4) go thru AAs in atlas & assign (based on strcmp & maxDist) to AA -> aaName

% (c) Jiri, Jan18, Jul18, Aug22

if nargin < 5
    maxDist = 0;
end
if nargin < 4
    addLateralization = false;
end

%% anatomy assignments from neurologists (A. Kalina)
if strcmp(atlasName, 'LabelsAK')
    ch_ass = H_ch.neurologyLabel;
    ch_dist = 0;
end

%% Yeo7 assignments
if strcmp(atlasName, 'yeo7_name')
    if isfield(H_ch, 'ass_yeo7_name')
        ch_ass = H_ch.ass_yeo7_name;
        ch_dist = H_ch.ass_yeo7_dist;
    end
end

%% Yeo17 assignments
if strcmp(atlasName, 'yeo17_name')
    if isfield(H_ch, 'ass_yeo17_name')
        ch_ass = H_ch.ass_yeo17_name;
        ch_dist = H_ch.ass_yeo17_dist;
    end
end

%% Destrieux atlas (no lateralization) assignments
if strcmp(atlasName, 'destrie_name')
    if isfield(H_ch, 'ass_destrie_name')
        ch_ass = H_ch.ass_destrie_name;
        ch_dist = H_ch.ass_destrie_dist;
    end
end

%% Destrieux atlas (wi lateralization) assignments
if strcmp(atlasName, 'destLat_name')
    if isfield(H_ch, 'ass_destLat_name')
        ch_ass = H_ch.ass_destLat_name;
        ch_dist = H_ch.ass_destLat_dist;
    end
end

%% Mars atlas (no lateralization) assignments
if strcmp(atlasName, 'mars_name')
    if isfield(H_ch, 'ass_mars_name')
        ch_ass = H_ch.ass_mars_name;
        ch_dist = H_ch.ass_mars_dist;
    end
end

%% TBD? cyto-architectonic assignments (motor plus sensory area combined)
if strcmp(atlasName, 'cytoarch_large_motorPlusSensory')
    ch_ass = H_ch.ass_cytoarchMap;
    ch_dist = 0;
end

%% TBD? cyto-architectonic assignments
if strcmp(atlasName, 'cytoarch_large')
    ch_ass = H_ch.ass_cytoarchMap;
    ch_dist = 0;
end

%% TO DO: cyto-architectonic assignments of individual areas in Anatomy toolbox

%% (2) in case of "multiple area" assignments (typically from neurologistLabels): 
%   - separate (separator = '/') to cell arrays 
%   - sort alphabetically
%   - remove 'wm'
%   - sort back
if ~isempty(strfind(ch_ass, '/'))
    ch_ass = str_splitSortDelWM(ch_ass);
end

%% (3) add laterality (L/R)? (based on MNI)
% !!! if x_MNI = 0 => assigned to R !!!
if addLateralization
    assert(exist('ch_ass','var') == 1);
    
    % laterality based on x MNI coordinates: L < 0 & R > 0
    if H_ch.MNI_x >= 0
        ch_side = 'R';
    elseif H_ch.MNI_x < 0
        ch_side = 'L';
    end
    ch_ass = [ch_side ' ' ch_ass];  % includes also empty space!
end

%% return the name of the anatomic area in the AA_list (first column)
aaName = 'n.a.';    % default
was_assigned = false;
for aa = 1:size(AA_list,1)
    for c = 1:size(AA_list{aa,2},2)
        if strcmpi(ch_ass, AA_list{aa,2}{c}) && ch_dist <= maxDist   % not case sensitive !!!
            if ~was_assigned
                aaName = AA_list{aa,1};
                was_assigned = true;
            else
                disp(['WARNING! channel already assigned to: ' aaName ', chAss = ' ch_ass ', name = ' H_ch.name]);
                error('channel already assigned, check the anatomy areas input table!');
            end
        end
    end
end
if ~was_assigned
    disp(['WARNING! channel not grouped: name = ' H_ch.name ', brainAtlas = ' atlasName ', ch-assignment = ' ch_ass ' (d=' num2str(ch_dist,2) ' mm)']);
end

