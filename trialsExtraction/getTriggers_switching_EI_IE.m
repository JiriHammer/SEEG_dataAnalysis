function getTriggers_switching_EI_IE(params)
% determines trials where there was a switch from internal to external task (and vice versa)
% determines trials which contained the selected keyword (in params.thisCategory)
% MUST match trials 'tr' with the same keyword in trials(tr).stim_cond
% saves result into respective column of 'trials_events' in each sessionCacheFile
%   - rows = trials
%   - values = zeros (event not found) and ones (event found)

% (c) Jiri, Mar16

display([' - finding events for: ' params.thisCategory]);

%% per session: load trials, find events and save
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load trials
    clear trials pred_raw
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    nTr = size(trials,1);
    trial_event = zeros(nTr, 1);  % default: no events found
    
    % find events
    for tr = 1:nTr-1
        if strcmp(params.thisCategory, 'E-I')
            if strcmp(trials(tr).stim_cond, 'E-task') && strcmp(trials(tr+1).stim_cond, 'I-task')
                trial_event(tr) = 1;
            end
        elseif strcmp(params.thisCategory, 'I-E')
            if strcmp(trials(tr).stim_cond, 'I-task') && strcmp(trials(tr+1).stim_cond, 'E-task')
                trial_event(tr) = 1;
            end
        elseif strcmp(params.thisCategory, 'E-E')
            if strcmp(trials(tr).stim_cond, 'E-task') && strcmp(trials(tr+1).stim_cond, 'E-task')
                trial_event(tr) = 1;
            end        
        elseif strcmp(params.thisCategory, 'I-I')
            if strcmp(trials(tr).stim_cond, 'I-task') && strcmp(trials(tr+1).stim_cond, 'I-task')
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
