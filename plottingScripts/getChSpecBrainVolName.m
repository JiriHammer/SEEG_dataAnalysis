function volumeName = getChSpecBrainVolName(channelStruct, volume2plot, groupInfo)
% TBD !!! (now = Nov22)
% returns volumeName (e.g. wDefault, wDorsal Attention, ...)
% based on channel information in channelStruct and specification in volume2plot

% (c) Jiri, Nov21

brainAtlas_name = strrep(volume2plot, 'ch-spec_', '');
volumeName = ch2roi_assignAA(brainAtlas_name, channelStruct, groupInfo);

%% old (now = Nov22): only for Yeo7 network
% volumeName = [];        % default
% if strcmp(volume2plot, 'ch-spec_Yeo7')
%     all_yeo7 = {'Visual';'Somatomotor';'Dorsal Attention';'Ventral Attention';'Limbic';'Frontoparietal';'Default'};
%     ass_yeo7 = channelStruct.ass_yeo7_name;
%     for nn = 1:size(all_yeo7,1)
%         if ~isempty(strfind(ass_yeo7, all_yeo7{nn}))
% %             volumeName = ['w' all_yeo7{nn}];    % old (added prefix 'w')
%             volumeName = all_yeo7{nn};
%         end
%     end
%     
%     if isempty(volumeName)
%         disp([volume2plot ': ch-specific brain volume not found: ass_yeo7 = ' ass_yeo7]);
%     end
% end
