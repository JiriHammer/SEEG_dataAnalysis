function getTriggers_ssep_LR(params)
% determines & saves trigger times and class labels

% (c) Jiri, May16

disp('triggering on left / right hand ...');

%% per session: load, get trigger times and assign labels, save
for sess = 1:size(params.storage.analyzedSess,2)
    
    t_trigger = [];
    clz_label = [];
    clz_names = {'LH', 1; 
                 'RH', 2};
    assert(ismember(params.thisCategory, clz_names(:,1)));  % must be either LH or RH
    
    % load trials = pred_raw
    clear trials
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    clear pred_raw;
    
    % get trigger times and assign labels
    for tr = 1:size(trials,1)
        
        % trigger times
        if strcmp(trials(tr).stim_type, params.thisCategory)
            t_trigger = cat(2, t_trigger, trials(tr).t_stimulus); 
            this_label = clz_names{ismember(clz_names(:,1), trials(tr).stim_type),2};   % 1 = LH, 2 = RH
            clz_label = cat(2, clz_label, this_label);
        end
        
%         % trigger times
%         if strcmp(trials(tr).stim_type, 'LH')
%             t_trigger = cat(2, t_trigger, trials(tr).t_stimulus); 
%             clz_label = cat(2, clz_label, 1);
%         end
    end
    
%     % RIGHT HAND: get trigger times and assign labels
%     for tr = 1:size(trials,1)
%         
%         % trigger times
%         if strcmp(trials(tr).stim_type, 'RH')
%             t_trigger = cat(2, t_trigger, trials(tr).t_stimulus); 
%             clz_label = cat(2, clz_label, 2);
%         end
%     end
    
    assert(size(t_trigger,2) == size(clz_label,2));
    
    % triggers
    clear triggers
    load(params.storage.sessionCacheFiles{sess}, 'triggers');   
    triggers = cat(2, triggers, cat(1, t_trigger, clz_label));
      
    % class names & labels
    classNamesLabels = clz_names;             % class names & labels
    
    % save
    save(params.storage.sessionCacheFiles{sess}, 'triggers','clz_names','classNamesLabels', '-append');
end

    