function esmType = getChannelType_ESM(chStruct_H)
% returns string of catted ESM categories (e.g. 'motor') based on ESM label
% ESM = electrical stimulation mapping
% category = higher order of abstraction that ESM label
% categories may overlap (-> esmType is catted into 1 string)
% chStruct_H = channel structure in H.channels(ch)
%    - should contain the field 'esm'
% ESM labels provided by neurologists
% output example: esmType = ' / handArmMotor / generalMotor' for ESM label 'L hand motor'

% (c) Jiri, Sep17, Sep18

esmType = [];

%% if there is no field 'esm' (el.stim. mapping)
if ~isfield(chStruct_H, 'esm')
    display(['WARNING: no ESM label for channel: ' chStruct_H.name ', signal type = ' chStruct_H.signalType]);
    esmType = [esmType 'n.a.'];
    return;
end
    
%% "hand/arm motor" channel group
esmLabels = {'hand motor', 'arm motor', 'shoulder motor'};
for lbl = 1:size(esmLabels,2)
    if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))
        esmType = [esmType ' / handArmMotor'];
    end
end
   
%% "eyes motor" channel group
esmLabels = {'eyes motor'};
for lbl = 1:size(esmLabels,2)
    if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))
        esmType = [esmType ' / eyesMotor'];
    end
end

%% "head motor" channel group
esmLabels = {'mouth motor', 'face motor', 'tongue motor', 'head motor'};
for lbl = 1:size(esmLabels,2)
    if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))
        esmType = [esmType ' / headMotor'];
    end
end

%% "motor" channel group (includes all motor responses: mouth, eyes, hands, legs, ...)
esmLabels = {'motor'};
for lbl = 1:size(esmLabels,2)
    if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))
        esmType = [esmType ' / generalMotor'];
    end
end

%% "sensory" channel group (includes all motor responses: mouth, eyes, hands, legs, ...)
esmLabels = {'sense', 'sensory'};
for lbl = 1:size(esmLabels,2)
    if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))
        esmType = [esmType ' / generalSensory'];
    end
end

%% "noResponse" channel group (used as "control group" until 27.09.2018, then renamed to "noResponse")
esmLabels = {'n.a.','n.s.r.'};
for lbl = 1:size(esmLabels,2)
    if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))
        %esmType = 'control';        % used as "control" until 27.09.2018
        esmType = [esmType ' / noResponse'];
    end
end




%% old "control group" (non-motor cortex) channels with areas exclusion
% commented out on 09.08.2018
% esmLabels = {'n.a.', 'n.s.r.'};
% aaAreas_exclude = {'Broca','MotorSensory','SPL/IPS'};
% for lbl = 1:size(esmLabels,2)
%     if ~isempty(strfind(chStruct_H.esm,esmLabels{lbl}))     % no (specific) ESM label
%         % define AA: 'thisAA' (based on the Anatomy toolbox)
%         elAss = se_TabList_jiri(chStruct_H.MNI_x, chStruct_H.MNI_y, chStruct_H.MNI_z);
%         AA_H_struct = chStruct_H.ass_cytoarchMap;
%         AA_toolbox = elAss.probabAss_area{1}; 
%         if ~strcmp(AA_H_struct, elAss.probabAss_area{1})
%             if strcmp(AA_toolbox, 'n.a.')
%                 thisAA = getLargerAnatomicalArea(AA_H_struct);
%             else
%                 thisAA = getLargerAnatomicalArea(AA_toolbox);
%             end
%         else
%             thisAA = getLargerAnatomicalArea(AA_toolbox);
%         end 
% 
%         % exclude areas (near) motor + sensory cortex channels
%         if ~ismember(thisAA, aaAreas_exclude)
%             esmType = 'control';
%         end
%     end
% end
