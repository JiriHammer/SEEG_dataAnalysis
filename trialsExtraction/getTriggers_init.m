function getTriggers_init(params, flag)
% initializes/allocates trigger labels to each session cache file
% based on the flag
% these variables are needed later to select the relevant trials

% (c) Jiri, Nov16

for sess = 1:size(params.storage.sessionCacheFiles,2)

    if strcmp(flag,'trials_events')
        clear pred_raw
        load(params.storage.sessionCacheFiles{sess}, 'pred_raw');   % load trials (pred_raw)
        nTr = size(pred_raw,1);        
        nCtg = params.nCtg;   % number of categories (incl. AND/OR)
        trials_events = nan(nTr, nCtg);     % trials x event categories
        save(params.storage.sessionCacheFiles{sess}, 'trials_events', '-append');
        
    elseif strcmp(flag,'trials_class')
        clear trials_events
        load(params.storage.sessionCacheFiles{sess}, 'trials_events');
        triggers_clz = trials_events(:,1);         % events of each class (for further merging by logical operators)
        save(params.storage.sessionCacheFiles{sess}, 'triggers_clz', '-append');
        
    elseif strcmp(flag,'times_labels')
        triggers = [];
        classNamesLabels = [];
        save(params.storage.sessionCacheFiles{sess}, 'triggers','classNamesLabels', '-append');
        
    else
        error(['uknown option: ' flag]);
    end

end
