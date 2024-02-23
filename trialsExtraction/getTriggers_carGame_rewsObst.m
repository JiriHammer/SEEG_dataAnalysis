function getTriggers_carGame_rewsObst(params)
% determines & saves trigger times and class labels

% (c) Jiri, May16

display('triggering on rewards / obstacles ...');

%% per session: load, get trigger times and assign labels, save
for sess = 1:size(params.storage.analyzedSess,2)
    
    t_trigger = [];
    clz_label = [];
    clz_names = {'rews', 1; 'obst', 2};
    
    % load
    clear trials
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    
    % REWARDS: get trigger times and assign labels
    for tr = 1:size(trials.rewards,1)
        
        % trigger times
        if trials.rewards(tr).wasCollected
            t_trigger = cat(2, t_trigger, trials.rewards(tr).timeCollected); 
            clz_label = cat(2, clz_label, 1);
        end
    end
    
    % OBSTACLES: get trigger times and assign labels
    for tr = 1:size(trials.penalty,1)
        
        % trigger times
        if trials.penalty(tr).wasCollected
            t_trigger = cat(2, t_trigger, trials.penalty(tr).timeCollected); 
            clz_label = cat(2, clz_label, 2);
        end
    end
    
    assert(size(t_trigger,2) == size(clz_label,2));
    triggers = cat(1, t_trigger, clz_label);
    
    % save
    save(params.storage.sessionCacheFiles{sess}, 'triggers','clz_names', '-append');
end

    