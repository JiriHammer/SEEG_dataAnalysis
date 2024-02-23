function getTriggers_gameType(params)
% extracts only defined paradigm from each session
% based on the keyword 'gameType' which should be loaded into sessionCache
% via params.load_D

% (c) Jiri, Nov16

display([' - finding events for: ' params.thisCategory]);

%% per session: load trials, find events and save
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load trials
    clear trials
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    nTr = size(trials,1);
    trial_event = zeros(nTr, 1);  % default: no events found
    
    % load paradigm type (=gameType)
    clear gameType
    load(params.storage.sessionCacheFiles{sess}, 'gameType');
    assert(exist('gameType','var') == 1);
    
    % find events
    for tr = 1:nTr
        if strcmp(params.thisCategory, gameType)
            trial_event(tr) = 1;
        end
    end
    
    % update trials-categories
    clear trials_events
    load(params.storage.sessionCacheFiles{sess}, 'trials_events');
    trials_events(:,params.thisCtg) = trial_event;
    save(params.storage.sessionCacheFiles{sess}, 'trials_events', '-append');
end
