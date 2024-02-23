function getTriggers_DES_intensity(params)
% determines trials which contained the selected keyword (in params.thisCategory)
% MUST match trials 'tr' with the same keyword in trials(tr).stim_type
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
    for tr = 1:nTr
        if strcmp(params.thisCategory, trials(tr).desIntensity)
            trial_event(tr) = 1;
        end
    end    
    
    % update trials-categories
    clear trials_events
    load(params.storage.sessionCacheFiles{sess}, 'trials_events');
    trials_events(:,params.thisCtg) = trial_event;
    save(params.storage.sessionCacheFiles{sess}, 'trials_events', '-append');
end
