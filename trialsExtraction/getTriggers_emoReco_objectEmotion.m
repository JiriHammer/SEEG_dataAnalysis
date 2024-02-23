function getTriggers_emoReco_objectEmotion(params)
% determines & saves trigger times and class labels

% (c) Jiri, Mar16

display([' - finding events for: ' params.thisCategory]);

%% per session: load trials, find events and save
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load trials
    clear trials
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    nTr = size(trials,1);
    trial_event = zeros(nTr, 1);  % default: no events found
    
    % find events
    for tr = 1:nTr
        if isfield(trials(tr).stimulus(1), 'obj_emo')
            if strcmp(params.thisCategory, trials(tr).stimulus(1).obj_emo)
                trial_event(tr) = 1;
            end
        end
    end
    
    % update trials-categories
    clear trials_events
    load(params.storage.sessionCacheFiles{sess}, 'trials_events');
    trials_events(:,params.thisCtg) = trial_event;
    save(params.storage.sessionCacheFiles{sess}, 'trials_events', '-append');
end

















%% OLD (22.11.2016): per session: load, get trigger times and assign labels, save
% for sess = 1:size(params.storage.analyzedSess,2)
%     
%     t_trigger = [];
%     clz_label = [];
%     clz_names = {
%         'emotions', 1; 
%         'objects', 2;
%         };
%     
%     % load
%     clear trials
%     load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
%     trials = pred_raw;
%     
%     % get trigger times and assign labels
%     for tr = 1:size(trials,1)
%         
%         % trigger times
%         if strcmp(params.triggering.cutPoint, 'onShowStimuli')
%             t_trigger = cat(2, t_trigger, trials(tr).stimulus(1).t_stim);   % hack: time of sound assumed same as those of objects!
%         elseif strcmp(params.triggering.cutPoint, 'onGetReady')
%             t_trigger = cat(2, t_trigger, trials(tr).t_ready);
%         elseif strcmp(params.triggering.cutPoint, 'onPause')
%             t_trigger = cat(2, t_trigger, trials(tr).t_pause);
%         elseif strcmp(params.triggering.cutPoint, 'onPromptAnswer')
%             t_trigger = cat(2, t_trigger, trials(tr).answer.t_choice);
%         elseif strcmp(params.triggering.cutPoint, 'onButtonPress')
%             t_trigger = cat(2, t_trigger, trials(tr).answer.t_answer);            
%         end
%                     
%         % class labels
%         if strcmp(trials(tr).stimulus(1).obj_emo, 'Emotions')
%             clz_label = cat(2, clz_label, 1);
%         elseif strcmp(trials(tr).stimulus(1).obj_emo, 'Objects')
%             clz_label = cat(2, clz_label, 2);
%         else
%             error(['uknown class label: ' trials(tr).stimulus(1).same_diff]);
%         end
%     end
%     assert(size(t_trigger,2) == size(clz_label,2));
%     triggers = cat(1, t_trigger, clz_label);
%     
%     % save
%     save(params.storage.sessionCacheFiles{sess}, 'triggers','clz_names', '-append');
% end
% 
%     