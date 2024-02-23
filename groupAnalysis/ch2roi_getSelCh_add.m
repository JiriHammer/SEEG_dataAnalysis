function addCh = ch2roi_getSelCh_add(i_ch_trials, ch_inds, chGroups, nn, cName)
% whether to add this channel to the selected ones or not
% !!! works only for 2 NN !!!
% 1) check if 'i_ch_trials' NOT already among the ch_inds
% 2) for NN = 2 (works only for 2 NN !!!)
%   3) in case a channel was assigned to the 1st group already 
%       (may happen if 2 different atlases are used: Yeo7 & Yeo17), then do not include!!!, works only for 2 NN !!!    
%   4) in case of BIP, check if one of the bipolar contacts is not in the other group 
%       (bip = B8-B9 & group = {'A1-A2', ..., 'B7-B8'}, contact B8 would contribute to both groups)

% (c) Jiri, Nov23

%% ---- version without checking of one contact from bipolar pair is also in the other group (good DTF result, Oct-Nov23)
% % init
% addCh = false;
% 
% % if NOT already among the ch_inds
% if ~ismember(i_ch_trials, ch_inds)
%     % in case a channel was assigned to the 1st group already, then do not include!!!, works only for 2 NN !!!
%     if (nn==1) || (nn == 2 && ~ismember(i_ch_trials, chGroups.ch_inds_trials{1,1})) 
%         addCh = true;
%     end
% end

%% ---- version with checking of one contact from bipolar pair is also in the other group
% init
addCh = false;

% check if NOT already among the ch_inds
if ~ismember(i_ch_trials, ch_inds)
    if nn==1         
        addCh = true;
    else
        % in case a channel was assigned to the 1st group already (may happen if 2 different atlases are used: Yeo7 & Yeo17), then do not include!!!, works only for 2 NN !!!    
        if ~ismember(i_ch_trials, chGroups.ch_inds_trials{1,1})
            % check if one of the bipolar contacts is not in the other group (bip = B8-B9 & group = {'A1-A2', ..., 'B7-B8'}, contact B8 would contribute to both groups)
            cName2 = strsplit(cName, '-');
            if size(cName2,2) == 1  % for CAR (no need to check the names)
                addCh = true;   
            else                    % for BIP
                g1_cNames = {};
                for k = 1:size(chGroups.ch_names{1},2)
                    g1_cNames = cat(2, g1_cNames, strsplit(chGroups.ch_names{1}{k}, '-'));
                end
                g1_cNames = unique(g1_cNames);
                if ~ismember(cName2{1}, g1_cNames) && ~ismember(cName2{2}, g1_cNames)
                    addCh = true;
                else
                    disp(['WARNING: ch = ' cName ', one of bipolar pair in the 1st group! Rejecting.']);
                end
            end
        end
    end
end

        